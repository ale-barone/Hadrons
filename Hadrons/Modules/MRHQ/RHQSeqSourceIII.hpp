/*
 * RHQSeqSourceIII.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MRHQ_RHQSeqSourceIII_hpp_
#define Hadrons_MRHQ_RHQSeqSourceIII_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Sequential gamma source                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MRHQ)

class RHQSeqSourceIIIPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RHQSeqSourceIIIPar,
                                    std::string,    q,
                                    unsigned int,   t,
                                    std::string,    mom,
                                    unsigned int,   index,
                                    Gamma::Algebra, gamma5,
                                    std::string,    gauge); 
};

template <typename FImpl, typename GImpl>
class TRHQSeqSourceIII: public Module<RHQSeqSourceIIIPar>
{
public:
    //FERM_TYPE_ALIASES(FImpl,);
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TRHQSeqSourceIII(const std::string name);
    // destructor
    virtual ~TRHQSeqSourceIII(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void makeSource(PropagatorField &src, const PropagatorField &q);
private:
    bool        hasPhase_{false};
    std::string momphName_, tName_;
};

MODULE_REGISTER_TMP(RHQSeqSourceIII, ARG(TRHQSeqSourceIII<FIMPL, GIMPL>), MRHQ);

/******************************************************************************
 *                     TRHQSeqSourceIII implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
TRHQSeqSourceIII<FImpl, GImpl>::TRHQSeqSourceIII(const std::string name)
: Module<RHQSeqSourceIIIPar>(name)
, momphName_ (name + "_momph")
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQSeqSourceIII<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().gauge};
    
    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQSeqSourceIII<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQSeqSourceIII<FImpl, GImpl>::setup(void)
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
void TRHQSeqSourceIII<FImpl, GImpl>::makeSource(PropagatorField &src, 
                                               const PropagatorField &q)
{
    
    auto &ph  = envGet(LatticeComplex, momphName_);
    auto &t   = envGet(Lattice<iScalar<vInteger>>, tName_);

    ComplexD           i(0.0,1.0);
    std::vector<RealD> p;
    p = strToVec<RealD>(par().mom);

    if (!hasPhase_)
    {
        envGetTmp(LatticeComplex, coor);
        ph = Zero();

        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }   

        ph = exp((RealD)(2*M_PI)*i*ph);
        LatticeCoordinate(t, Tp);
        hasPhase_ = true;
    }
    
    auto &field = envGet(PropagatorField, par().q);
    const auto &gaugefield = envGet(GaugeField, par().gauge);
    const auto &index = par().index;
    const auto gauge_x = peekLorentz(gaugefield, 0);
    const auto gauge_y = peekLorentz(gaugefield, 1);
    const auto gauge_z = peekLorentz(gaugefield, 2);

    Gamma g5(par().gamma5);
    Gamma gx(Gamma::Algebra::GammaX);
    Gamma gy(Gamma::Algebra::GammaY);
    Gamma gz(Gamma::Algebra::GammaZ);

    const PropagatorField Dx = GImpl::CovShiftForward(gauge_x,0,field) - GImpl::CovShiftBackward(gauge_x,0,field);
    const PropagatorField Dy = GImpl::CovShiftForward(gauge_y,1,field) - GImpl::CovShiftBackward(gauge_y,1,field);
    const PropagatorField Dz = GImpl::CovShiftForward(gauge_z,2,field) - GImpl::CovShiftBackward(gauge_z,2,field);

    Gamma::Algebra gi; 
    switch(par().index){
        case 0:
            gi = Gamma::Algebra::GammaX;
            break;
        case 1:
            gi = Gamma::Algebra::GammaY;
            break;
        case 2:
            gi = Gamma::Algebra::GammaZ;
            break;
        case 3:
            gi = Gamma::Algebra::GammaT;
            break;
        default:
            HADRONS_ERROR(Argument, "Index must be in {0, 1, 2, 3}."); 
    }

    src = gi*gx*g5 * Dx - gx*gi*g5 * Dx
          + gi*gy*g5 * Dy - gy*gi*g5 * Dy
          + gi*gz*g5 * Dz - gz*gi*g5 * Dz;
        
    src = where((t == par().t), ph*src, 0.*src);
    src = 0.5*src;

}

template <typename FImpl, typename GImpl>
void TRHQSeqSourceIII<FImpl, GImpl>::execute(void)
{

    LOG(Message) << "Generating "
                 << "sequential source(s) at t=" << par().t << std::endl;


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

#endif // Hadrons_MRHQ_RHQSeqSourceIII_hpp_
