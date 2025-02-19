#pragma once
#ifndef LI_LeptonDepthFunction_H
#define LI_LeptonDepthFunction_H

#include <memory>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/distributions/primary/vertex/DepthFunction.h"

namespace LI {
namespace dataclassses {
class InteractionSignature;
}
namespace distributions {

class LeptonDepthFunction : virtual public DepthFunction {
friend cereal::access;
private:
    double mu_alpha = 1.76666667e-3;
    double mu_beta = 2.0916666667e-6;
    double tau_alpha = 1.473684210526e1;
    double tau_beta = 2.6315789473684212e-7;
    double scale = 1.0;
    double max_depth = 3e7;
    std::set<LI::dataclasses::Particle::ParticleType> tau_primaries = {LI::dataclasses::Particle::ParticleType::NuTau, LI::dataclasses::Particle::ParticleType::NuTauBar};
public:
    LeptonDepthFunction();
    double GetLeptonDepthFunctionReturnValue(LI::dataclasses::InteractionSignature const & signature, double energy) const; 
    void SetMuParams(double mu_alpha, double mu_beta);
    void SetTauParams(double tau_alpha, double tau_beta);
    void SetScale(double scale);
    void SetMaxDepth(double max_depth);
    double GetMuAlpha() const;
    double GetMuBeta() const;
    double GetTauAlpha() const;
    double GetTauBeta() const;
    double GetScale() const;
    double GetMaxDepth() const;
    virtual double operator()(LI::dataclasses::InteractionSignature const & signature, double energy) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
        } else {
            archive(::cereal::make_nvp("MuAlpha", mu_alpha));
            archive(::cereal::make_nvp("MuBeta", mu_beta));
            archive(::cereal::make_nvp("TauAlpha", tau_alpha));
            archive(::cereal::make_nvp("TauBeta", tau_beta));
            archive(::cereal::make_nvp("Scale", scale));
            archive(::cereal::make_nvp("MaxDepth", max_depth));
            archive(::cereal::make_nvp("TauPrimaries", tau_primaries));
            throw std::runtime_error("LeptonDepthFunction only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            archive(::cereal::make_nvp("MuAlpha", mu_alpha));
            archive(::cereal::make_nvp("MuBeta", mu_beta));
            archive(::cereal::make_nvp("TauAlpha", tau_alpha));
            archive(::cereal::make_nvp("TauBeta", tau_beta));
            archive(::cereal::make_nvp("Scale", scale));
            archive(::cereal::make_nvp("MaxDepth", max_depth));
            archive(::cereal::make_nvp("TauPrimaries", tau_primaries));
            throw std::runtime_error("LeptonDepthFunction only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(DepthFunction const & distribution) const override;
    virtual bool less(DepthFunction const & distribution) const override;
};


} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::LeptonDepthFunction, 0);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::DepthFunction, LI::distributions::LeptonDepthFunction);

#endif // LI_LeptonDepthFunction_H
