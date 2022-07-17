/*
 * RHQInsertionIV.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MRHQ_RHQInsertionIV_hpp_
#define Hadrons_MRHQ_RHQInsertionIV_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            RHQInsertionIV                                  *
 ******************************************************************************/
GRID_SERIALIZABLE_ENUM(OpIVFlag, undef, Chroma, 0, LeftRight, 1);

BEGIN_MODULE_NAMESPACE(MRHQ)

class RHQInsertionIVPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RHQInsertionIVPar,
                                    std::string,    q,
                                    unsigned int,   index,
                                    Gamma::Algebra, gamma5,
                                    std::string,    gauge,
                                    std::string,    propTwist,
                                    OpIVFlag,       flag);
};

template <typename FImpl, typename GImpl>
class TRHQInsertionIV: public Module<RHQInsertionIVPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TRHQInsertionIV(const std::string name);
    // destructor
    virtual ~TRHQInsertionIV(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RHQInsertionIV, ARG(TRHQInsertionIV<FIMPL, GIMPL>), MRHQ);

/******************************************************************************
 *                       TRHQInsertionIV implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
TRHQInsertionIV<FImpl, GImpl>::TRHQInsertionIV(const std::string name)
: Module<RHQInsertionIVPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionIV<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().gauge};
    
    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionIV<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionIV<FImpl, GImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionIV<FImpl, GImpl>::execute(void)
{

    LOG(Message) << "Applying Improvement term IV with index " << par().index
                 << " and gamma5=" << par().gamma5 
                 << " to '" << par().q 
                 << "' with flag '" << par().flag << "'"
                 << std::endl;
    
    if (par().gamma5 != Gamma::Algebra::Gamma5 && par().gamma5 != Gamma::Algebra::Identity)
    {
        HADRONS_ERROR(Argument, "gamma5 must be either 'Gamma5' or 'Identity'."); 
    }
    Gamma g5(par().gamma5);

    // In case the propagator we are applying the RHQInseriton on contains twist
    ComplexD          i(0.0,1.0);

    ComplexD ph_f_x = exp(-i*0.);
    ComplexD ph_b_x = exp(i*0.);
    
    ComplexD ph_f_y = exp(-i*0.);
    ComplexD ph_b_y = exp(i*0.);

    ComplexD ph_f_z = exp(-i*0.);
    ComplexD ph_b_z = exp(i*0.);

    if (!par().propTwist.empty())
    {
        std::vector<RealD> p;
        p = strToVec<RealD>(par().propTwist);
    
        double q_x    = (2*M_PI*p[0]/env().getDim(0));
        double q_y    = (2*M_PI*p[1]/env().getDim(1));
        double q_z    = (2*M_PI*p[2]/env().getDim(2));
        
        ph_f_x = exp(i*q_x);
        ph_b_x = exp(-i*q_x);

        ph_f_y = exp(i*q_y);
        ph_b_y = exp(-i*q_y);

        ph_f_z = exp(i*q_z);
        ph_b_z = exp(-i*q_z);
    }

    
    auto &field = envGet(PropagatorField, par().q);
    const auto &gaugefield = envGet(GaugeField, par().gauge);
    const auto gauge_x = peekLorentz(gaugefield, 0);
    const auto gauge_y = peekLorentz(gaugefield, 1);
    const auto gauge_z = peekLorentz(gaugefield, 2);

    Gamma gx(Gamma::Algebra::GammaX);
    Gamma gy(Gamma::Algebra::GammaY);
    Gamma gz(Gamma::Algebra::GammaZ);

    const PropagatorField Dx = ph_f_x*GImpl::CovShiftForward(gauge_x,0,field) - ph_b_x*GImpl::CovShiftBackward(gauge_x,0,field);
    const PropagatorField Dy = ph_f_y*GImpl::CovShiftForward(gauge_y,1,field) - ph_b_y*GImpl::CovShiftBackward(gauge_y,1,field);
    const PropagatorField Dz = ph_f_z*GImpl::CovShiftForward(gauge_z,2,field) - ph_b_z*GImpl::CovShiftBackward(gauge_z,2,field);

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
    
    auto &out = envGet(PropagatorField, getName());
    if (par().flag == OpIVFlag::Chroma)
    {     
        PropagatorField insertion =
            gx*g5*gi * Dx 
          + gy*g5*gi * Dy
          + gz*g5*gi * Dz;
        
        out = insertion;
    }
    else if (par().flag == OpIVFlag::LeftRight)
    {        
        PropagatorField insertion = 
            gi*gx*g5 * Dx - gx*gi*g5 * Dx
          + gi*gy*g5 * Dy - gy*gi*g5 * Dy
          + gi*gz*g5 * Dz - gz*gi*g5 * Dz;

        out = -0.5*insertion;
    }
 
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MRHQ_RHQInsertionIV_hpp_
