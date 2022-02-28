/*
 * RHQSeqSourceIV.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MRHQ_RHQSeqSourceIV_hpp_
#define Hadrons_MRHQ_RHQSeqSourceIV_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Sequential gamma source                            *
 ******************************************************************************/
GRID_SERIALIZABLE_ENUM(OpIVMomType, undef, Sink, 0, Twist, 1);

BEGIN_MODULE_NAMESPACE(MRHQ)

class RHQSeqSourceIVPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RHQSeqSourceIVPar,
                                    std::string,    q,
                                    unsigned int,   t,
                                    std::string,    mom,
                                    unsigned int,   index,
                                    Gamma::Algebra, gamma5,
                                    std::string,    gauge,
                                    OpIVMomType,    momType); 
};

template <typename FImpl, typename GImpl>
class TRHQSeqSourceIV: public Module<RHQSeqSourceIVPar>
{
public:
    //FERM_TYPE_ALIASES(FImpl,);
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TRHQSeqSourceIV(const std::string name);
    // destructor
    virtual ~TRHQSeqSourceIV(void) {};
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

MODULE_REGISTER_TMP(RHQSeqSourceIV, ARG(TRHQSeqSourceIV<FIMPL, GIMPL>), MRHQ);

/******************************************************************************
 *                     TRHQSeqSourceIV implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
TRHQSeqSourceIV<FImpl, GImpl>::TRHQSeqSourceIV(const std::string name)
: Module<RHQSeqSourceIVPar>(name)
, momphName_ (name + "_momph")
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQSeqSourceIV<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().gauge};
    
    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQSeqSourceIV<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQSeqSourceIV<FImpl, GImpl>::setup(void)
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
void TRHQSeqSourceIV<FImpl, GImpl>::makeSource(PropagatorField &src, 
                                               const PropagatorField &q)
{
    
    auto &ph  = envGet(LatticeComplex, momphName_);
    auto &t   = envGet(Lattice<iScalar<vInteger>>, tName_);

    Complex           i(0.0,1.0);
    std::vector<Real> p;
    p = strToVec<Real>(par().mom);

    if (!hasPhase_)
    {
        envGetTmp(LatticeComplex, coor);
        ph = Zero();

        if (par().momType == OpIVMomType::Sink)
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
    const auto gauge_x = peekLorentz(gaugefield, 0);
    const auto gauge_y = peekLorentz(gaugefield, 1);
    const auto gauge_z = peekLorentz(gaugefield, 2);

    Gamma g5(par().gamma5);
    Gamma gx(Gamma::Algebra::GammaX);
    Gamma gy(Gamma::Algebra::GammaY);
    Gamma gz(Gamma::Algebra::GammaZ);

    const PropagatorField Dx_f = GImpl::CovShiftForward(gauge_x,0,field);
    const PropagatorField Dy_f = GImpl::CovShiftForward(gauge_y,1,field);
    const PropagatorField Dz_f = GImpl::CovShiftForward(gauge_z,2,field);
    const PropagatorField Dx_b = GImpl::CovShiftBackward(gauge_x,0,field);
    const PropagatorField Dy_b = GImpl::CovShiftBackward(gauge_y,1,field);
    const PropagatorField Dz_b = GImpl::CovShiftBackward(gauge_z,2,field);

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

    PropagatorField src_f = 
        gi*gx*g5 * Dx_f - gx*gi*g5 * Dx_f
      + gi*gy*g5 * Dy_f - gy*gi*g5 * Dy_f
      + gi*gz*g5 * Dz_f - gz*gi*g5 * Dz_f;

    PropagatorField src_b =
        gi*gx*g5 * Dx_b - gx*gi*g5 * Dx_b
      + gi*gy*g5 * Dy_b - gy*gi*g5 * Dy_b
      + gi*gz*g5 * Dz_b - gz*gi*g5 * Dz_b;

    if (index == 3)
    {        
        // for forward shift
        int t_f = par().t-1;
        src_f = where((t == t_f), ph*src_f, 0.*src_f);
        // for backward shift
        int t_b = par().t+1;
        src_b = where((t == t_b), ph*src_b, 0.*src_b);

	    src = src_b - src_f;
    }
    else
    {        
        double ph_n = (2*M_PI*p[index]/env().getDim(index));
		ComplexD ph_f = exp(i*ph_n);
		ComplexD ph_b = exp(-i*ph_n);

	    src = -ph_f*src_f + ph_b*src_b;
        src = where((t == par().t), ph*src, 0.*src);
    }

}

template <typename FImpl, typename GImpl>
void TRHQSeqSourceIV<FImpl, GImpl>::execute(void)
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

#endif // Hadrons_MRHQ_RHQSeqSourceIV_hpp_
