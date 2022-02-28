/*
 * RHQSeqSourceII.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
 * Author: Alessandro Barone <barone1618@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */
 
/*  END LEGAL */

#ifndef Hadrons_MRHQ_RHQSeqSourceII_hpp_
#define Hadrons_MRHQ_RHQSeqSourceII_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source
 -----------------------------
 * src_x = q_x * theta(x_3 - tA) * theta(tB - x_3) * gamma * exp(i x.mom)
 
 * options:
 - q: input propagator (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - gamma: gamma product to insert (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         Sequential gamma source                            *
 ******************************************************************************/
GRID_SERIALIZABLE_ENUM(OpIIMomType, undef, Sink, 0, Twist, 1);

BEGIN_MODULE_NAMESPACE(MRHQ)

class RHQSeqSourceIIPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RHQSeqSourceIIPar,
                                    std::string,    q,
                                    unsigned int,   t, 
                                    std::string,    mom,
                                    OpIIMomType,    momType,
                                    std::string,    gauge,
                                    unsigned int,   index,
                                    Gamma::Algebra, gamma5); 
};

template <typename FImpl, typename GImpl>
class TRHQSeqSourceII: public Module<RHQSeqSourceIIPar>
{
public:
    //FERM_TYPE_ALIASES(FImpl,);
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TRHQSeqSourceII(const std::string name);
    // destructor
    virtual ~TRHQSeqSourceII(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void makeSourcePhase(PropagatorField &src, const PropagatorField &q);
    void makeSource(PropagatorField &src, const PropagatorField &q);
private:
    bool        hasPhase_{false};
    std::string momphName_, tName_;
};

MODULE_REGISTER_TMP(RHQSeqSourceII, ARG(TRHQSeqSourceII<FIMPL, GIMPL>), MRHQ);

/*******************************************************************************
 *                     TRHQSeqSourceII implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
TRHQSeqSourceII<FImpl, GImpl>::TRHQSeqSourceII(const std::string name)
: Module<RHQSeqSourceIIPar>(name)
, momphName_ (name + "_momph")
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQSeqSourceII<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().gauge};
    
    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQSeqSourceII<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQSeqSourceII<FImpl, GImpl>::setup(void)
{  
    if (envHasType(PropagatorField, par().q))
    {
        envCreateLat(PropagatorField, getName());
    }
    else
    {
        HADRONS_ERROR_REF(ObjectType, "object '" + par().q 
                          + "' has an incompatible type ("
                          + env().getObjectType(par().q)
                          + ")", env().getObjectAddress(par().q))
    }
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envCacheLat(LatticeComplex, momphName_);
    envTmpLat(LatticeComplex, "coor");
}

// execution ///////////////////////////////////////////////////////////////////
// makeSource //////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQSeqSourceII<FImpl, GImpl>::makeSource(PropagatorField &src, 
                                               const PropagatorField &q)
{
    
    auto &ph  = envGet(LatticeComplex, momphName_);
    auto &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    Gamma g5(par().gamma5);

    Complex           i(0.0,1.0);
    std::vector<Real> p;
    p = strToVec<Real>(par().mom);

    if (!hasPhase_)
    {
        envGetTmp(LatticeComplex, coor);
        ph = Zero();

        if (par().momType == OpIIMomType::Sink)
        {
            for(unsigned int mu = 0; mu < env().getNd(); mu++)
            {
                LatticeCoordinate(coor, mu);
                ph = ph + (p[mu]/env().getDim(mu))*coor;
            }   
        }

        ph = exp((Real)(2*M_PI)*i*ph);
        LatticeCoordinate(t, Tp);
        hasPhase_ = true;
    }

    auto &field = envGet(PropagatorField, par().q);
    const auto &gaugefield = envGet(GaugeField, par().gauge);
    const auto &index = par().index;
    const auto internal_gauge = peekLorentz(gaugefield, index);
    
    PropagatorField src_f = GImpl::CovShiftForward(internal_gauge, index, field);
    PropagatorField src_b = GImpl::CovShiftBackward(internal_gauge, index, field);

    if (index == 3)
    {        
        // for forward shift
        int t_f = par().t-1;
        src_f = where((t == t_f), ph*(g5*src_f), 0.*src_f);
        // for backward shift
        int t_b = par().t+1;
        src_b = where((t == t_b), ph*(g5*src_b), 0.*src_b);

	    src = src_b - src_f;
    }
    else
    {        
        double ph_n = (2*M_PI*p[index]/env().getDim(index));
		ComplexD ph_f = exp(i*ph_n);
		ComplexD ph_b = exp(-i*ph_n);

	    src = -ph_f*src_f + ph_b*src_b;
        src = where((t == par().t), ph*(g5*src), 0.*src);
    }

}

// execution//// ///////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQSeqSourceII<FImpl, GImpl>::execute(void)
{

    LOG(Message) << "Generating RHQ improved sequential source(s) "
                 << "with index" << par().index
                 << " at t=" << par().t << std::endl;

    if (envHasType(PropagatorField, par().q))
    {
        auto  &src = envGet(PropagatorField, getName()); 
        auto  &q   = envGet(PropagatorField, par().q);

        LOG(Message) << "Using propagator '" << par().q << "'" << std::endl;
        makeSource(src, q);
    }
    else
    {
        auto  &src = envGet(std::vector<PropagatorField>, getName()); 
        auto  &q   = envGet(std::vector<PropagatorField>, par().q);

        for (unsigned int i = 0; i < q.size(); ++i)
        {
            LOG(Message) << "Using element " << i << " of propagator vector '" 
                         << par().q << "'" << std::endl;
            makeSource(src[i], q[i]);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MRHQ_RHQSeqSourceII_hpp_
