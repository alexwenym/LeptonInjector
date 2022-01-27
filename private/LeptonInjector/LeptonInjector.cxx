#include <cassert>
#include <fstream>
#include <algorithm>

#include "LeptonInjector/LeptonInjector.h"
#include "LeptonInjector/EventProps.h"

#include "phys-services/LICrossSection.h"

#include "earthmodel-service/Path.h"

// namespace constants = boost::math::constants;

namespace LeptonInjector{


//--------------
//Config objects



//-----------
//Module base

LI_Position SampleFromDisk(std::shared_ptr<LI_random> rand, double radius, double zenith, double azimuth) {
    //choose a random point on a disk laying in the xy plane
    double t = rand->Uniform(0, 2 * Constants::pi);
    double u = rand->Uniform() + rand->Uniform();
    double r = ((u > 1.) ? (2. - u) : u) * radius;
    LI_Position  pos(r * cos(t) ,r * sin(t), 0.0);
    //now rotate to make the disc perpendicular to the requested normal vector
    pos = RotateY(pos, zenith);
    pos = RotateZ(pos, azimuth);
    return(pos);
}

double SampleEnergy(std::shared_ptr<LI_random> rand, double energyMinimum, double energyMaximum, double powerlawIndex) {
    if(energyMinimum == energyMaximum)
        return energyMinimum; //return the only allowed energy

    if(powerlawIndex == 1.0) //sample uniformly in log space
        return std::pow(10.0, rand->Uniform(log10(energyMinimum),log10(energyMaximum)));
    else{
        double u = rand->Uniform();
        double energyP = (1-u) * std::pow(energyMinimum, 1-powerlawIndex) + u * std::pow(energyMaximum, 1-powerlawIndex);
        return std::pow(energyP, 1/(1-powerlawIndex));
    }
}

// this function returns a pair of angles
std::pair<double,double> computeFinalStateAngles(Particle::ParticleType finalType1, Particle::ParticleType finalType2, double target_mass, double E_total, double x, double y) {
    double theta1=0, theta2=0;

    //first particle is a lepton, which covers CC, NC, and leptonic GR
    if(isLepton(finalType1)) {
        double m1=Particle(finalType1).GetMass();
        double E1 = (1 - y) * E_total;
        double cos_theta1, kE1;

        if(!isLepton(finalType2)){ //CC and NC have Hadrons as particle 2
            //squared kinetic energy of final state particle 1:
            double kE1sq=E1*E1 - m1*m1;
            if(kE1sq<=0){
                throw std::runtime_error("Negative kinetic energy. Not good");
            }
            cos_theta1=(E1 - x * y * target_mass - m1*m1/(2*E_total))/sqrt(kE1sq);
            kE1=sqrt(kE1sq);
        }
        else{ //leptonic GR
            double m_e = Constants::electronMass;

            if(E1<=0){ throw std::runtime_error("Bjorken Y > 1?"); }


            cos_theta1=1 - (m_e*m_e + 2*m_e*E_total - m1*m1)/(2*E_total*E1);
            kE1=E1;
        }

        if(cos_theta1<-1){
            // commented out until new logger is implemented
            //        log_warn_stream("cos(theta) underflow (" << cos_theta1 << "); rounding up to -1"
            //          "\n(E_total=" << E_total/I3Units::GeV << " x=" << x << " y=" << y << ")");
            cos_theta1=-1;
        }
        else if(cos_theta1>1){
            //tell the user if the difference was large enough to plausibly not be just round-off

            // need new logger
            //        if((cos_theta1-1)>1e-3)
            //          log_warn_stream("cos(theta) overflow (" << cos_theta1 << "); rounding down to 1"
            //            "\n(E_total=" << E_total/I3Units::GeV << " x=" << x << " y=" << y <<")");
            cos_theta1=1;
        }

        theta1=acos(cos_theta1);

        //longitudinal component of final state particle 2 momentum:
        double p_long=E_total-kE1*cos_theta1;
        //transverse component of final state particle 2 momentum:
        double p_trans=kE1*sin(theta1);
        theta2=atan(p_trans/p_long);
    }
    //otherwise we have hadronic GR, so both final state masses are unknown
    //and there isn't much we can do, so leave everything colinear
    // angle of particle 1 from initial dir
    return(std::make_pair(theta1,theta2));
}

// take a direction, deflect that direction by a distance /zenith/
//    rotate the new direction around the initial direction by /azimuth/
// So the zenith and azimuth are only what their names would suggest in the coordinate system where
//    /base/ is the \hat{z} axis

void FillTree(std::shared_ptr<LI_random> rand, Particle::ParticleType finalType1, Particle::ParticleType finalType2, LI_Position vertex, LI_Direction dir, double target_mass, LICrossSection::finalStateRecord const &fs, double energy, BasicEventProperties& properties, std::array<h5Particle,3>& particle_tree) {

    std::pair<double,double> relativeZeniths = computeFinalStateAngles(finalType1, finalType2, target_mass, energy, fs.x, fs.y);
    double azimuth1 = rand->Uniform(0, 2 * Constants::pi);
    double azimuth2 = azimuth1 + (azimuth1<Constants::pi ? 1 : -1)*Constants::pi ;

    (particle_tree)[0]=  h5Particle(true,
            static_cast<int32_t>(deduceInitialType(finalType1, finalType2)),
            vertex,
            dir,
            energy
            );

    //Make the first final state particle
    (particle_tree)[1] = h5Particle( false,
            static_cast<int32_t>(finalType1),
            vertex,
            rotateRelative(dir,relativeZeniths.first,azimuth1),
            kineticEnergy(finalType1,(1-fs.y)*energy)
            );

    (particle_tree)[2] = h5Particle(false,
            static_cast<int32_t>(finalType2),
            vertex,
            rotateRelative(dir,relativeZeniths.second,azimuth2),
            kineticEnergy(finalType2,fs.y*energy)
            );

    properties.totalEnergy=energy;
    properties.zenith=dir.zenith;
    properties.azimuth=dir.azimuth;
    properties.finalStateX=fs.x;
    properties.finalStateY=fs.y;
    properties.finalType1= static_cast<int32_t>(finalType1);
    properties.finalType2= static_cast<int32_t>(finalType2);
    properties.initialType=static_cast<int32_t>(deduceInitialType(finalType1, finalType2));
    properties.x = vertex.GetX();
    properties.y = vertex.GetY();
    properties.z = vertex.GetZ();
}

InteractionRecord RangedLeptonInjector::GenerateEvent() {
    InteractionRecord event = NewRecord();

    // Choose a target momentum
    target_momentum_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose an energy
    energy_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a direction on the sphere
    direction_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a position for the vertex
    position_distribution->Sample(random, earth_model, cross_sections, event);

    // Sample the cross section and final state
    SampleCrossSection(event);
    
    // Sample decay angle of photon
    SampleSecondaryDecay(event);
    
    // Sample pair production location 
    SamplePairProduction(event);
    
    injected_events += 1;
    return event;
}


//-----------------------
//Volume injection module

InteractionRecord VolumeLeptonInjector::GenerateEvent() {
    InteractionRecord event = NewRecord();

    // Choose a target momentum
    target_momentum_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose an energy
    energy_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a direction on the sphere
    direction_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a position for the vertex
    position_distribution->Sample(random, earth_model, cross_sections, event);

    // Sample the cross section and final state
    SampleCrossSection(event);
    injected_events += 1;
    return event;
}

} // namespace LeptonInjector

