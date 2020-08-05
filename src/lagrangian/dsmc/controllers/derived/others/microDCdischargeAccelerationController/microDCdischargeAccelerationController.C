/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "microDCdischargeAccelerationController.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "vectorList.H"
#include <math.h>
#include <unistd.h>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(microDCdischargeAccelerationController, 0);

    addToRunTimeSelectionTable(dsmcStateController, microDCdischargeAccelerationController, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    microDCdischargeAccelerationController::microDCdischargeAccelerationController
    (
        Time& t,
        dsmcCloud& cloud,
        const dictionary& dict
    )
    :
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
   
    // Pressure measurement
    nTimeSteps_(0.0),
    rhoNMean_(controlZone().size(), scalar(0.0)),
    rhoNMeanXnParticle_(controlZone().size(), scalar(0.0)),
    rhoMMeanXnParticle_(controlZone().size(), scalar(0.0)),
    momentumMeanXnParticle_(controlZone().size(), vector::zero),
    linearKEMeanXnParticle_(controlZone().size(), scalar(0.0)),
    rhoN_(controlZone().size(), scalar(0.0)),
    rhoNMean2_(controlZone().size(), scalar(0.0)),
    rhoMMean_(controlZone().size(), scalar(0.0)),
    linearKEMean_(controlZone().size(), scalar(0.0)),
    UMean_(controlZone().size(), vector::zero),
    translationalT_(controlZone().size(), scalar(0.0)),
    p_(controlZone().size(), scalar(0.0)),
    densityOnly_(1),
    controlVolume_(0.0),
    measuredParcels_(0.0),
    avParcelDensity_(0.0),
    pressuresum_(0.0),
    pressureavrg_(0.0),
    particlesToAccelerate_(mesh_.nCells()),
    accelerations_(mesh_.nCells(), vector::zero),
    
    //Cell ordering by (i,j) 
    
    cellzoneStart_(propsDict_.lookup("startPoint")),
    cellzoneEnd_(propsDict_.lookup("endPoint")),
    x_(readScalar(propsDict_.lookup("number_of_cells_in_xZone"))),   
    y_(readScalar(propsDict_.lookup("number_of_cells_in_yZone"))), 
       
    //Electric Source properties
    voltage_(readScalar(propsDict_.lookup("voltage"))),
    
    //Micro-gaps breakdown voltage properties
     workFunction_(readScalar(propsDict_.lookup("workFunction"))),
     enhancement_(readScalar(propsDict_.lookup("enhancementFactor"))),
     sEE_(readScalar(propsDict_.lookup("SEECoeff"))), 
     sEEm_(0.0),
     rhoI_(0.0),
     current_rhoI_(controlZone().size(), scalar(0.0)),
     rhoNMeanXnIonParticle_(controlZone().size(), scalar(0.0)),
     rhoNMeanXnIonParticle2_(controlZone().size(), scalar(0.0)),
     
     //Electrode's geometry
      gap_(readScalar(propsDict_.lookup("gap"))),
      ncathode_(readScalar(propsDict_.lookup("number_of_cathodes"))),
      nanode_(readScalar(propsDict_.lookup("number_of_anodes"))),
      cathode1xo_(propsDict_.lookup("cathode1_start_point")),
      cathode1xi_(propsDict_.lookup("cathode1_end_point")),
      cathode2xo_(propsDict_.lookup("cathode2_start_point")),
      cathode2xi_(propsDict_.lookup("cathode2_end_point")),
      anode1xo_(propsDict_.lookup("anode1_start_point")),
      anode1xi_(propsDict_.lookup("anode1_end_point")),
      anode2xo_(propsDict_.lookup("anode2_start_point")),
      anode2xi_(propsDict_.lookup("anode2_end_point")),   
     ionCounter_(0.0)
    {
 
        writeInTimeDir_ = false;
        writeInCase_ = false;
        singleValueController() = true;
        
        
         // standard to reading typeIds ------------ 
    const List<word> molecules (propsDict_.lookup("typeIds"));

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcInflowPatch::dsmcInflowPatch()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
        
    }
    
    Info << "controlZone().size() = " << controlZone().size() << endl;
 
  } 
    
    
    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

microDCdischargeAccelerationController::~microDCdischargeAccelerationController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void microDCdischargeAccelerationController::initialConfiguration()
{
}

void microDCdischargeAccelerationController::calculateProperties()
{
}  
   
   
void microDCdischargeAccelerationController::controlParcelsBeforeMove()
{ 

  //********Pressure measurement ***********************
    Info <<  " microDCdischargeAccelerationController successfully on!!! "  << nl << endl;  

    const labelList& cells = mesh_.cellZones()[regionId_];
    const List< DynamicList<dsmcParcel*> >& cellOccupancy= cloud_.cellOccupancy();
    nTimeSteps_ += 1.0;

    if(time_.samplingTime())
    {
        forAll(controlZone(), c)
        {
            const label& cellI = controlZone()[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
            
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];

                if(findIndex(typeIds_, p->typeId()) != -1)
                {        
                    const scalar& mass = cloud_.constProps(p->typeId()).mass();
                    const scalar& nParticle = cloud_.nParticle();
                    
//                     if(mass > 0)
//                     {
//                         mass_ = mass;
//                     }
                    
                    //rhoNMean_[cellI]++;
                     
                
                    if(cloud_.axisymmetric())
                    {
                            const point& cC = cloud_.mesh().cellCentres()[cellI];
                            
                            scalar radius = cC.y();
                            
                            scalar RWF = 1.0;
                            
                            RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
                        
                            rhoNMeanXnParticle_[cellI] += RWF*nParticle;
                            rhoMMeanXnParticle_[cellI] += mass*RWF*nParticle;
                            momentumMeanXnParticle_[cellI] += mass*p->U()*RWF*nParticle;
                            linearKEMeanXnParticle_[cellI] += mass*(p->U() & p->U())*RWF*nParticle;
                    }
                    else
                    {
                            rhoNMeanXnParticle_[cellI] += nParticle;
                            rhoMMeanXnParticle_[cellI] += mass*nParticle;
                            momentumMeanXnParticle_[cellI] += mass*p->U()*nParticle;
                            linearKEMeanXnParticle_[cellI] += mass*(p->U() & p->U())*nParticle;          
                    }
                } 
            }
        }
    } 
       
    if(time_.averagingTime())
    {
         
        const scalar& nAvTimeSteps = nTimeSteps_;

        forAll(rhoNMean_, c)
        {
            const label& cellI = controlZone()[c];
            const scalar& cellVolume = mesh_.cellVolumes()[cellI];
            
            if(rhoNMeanXnParticle_[c] > VSMALL)
            {
                rhoNMean_[c] =  rhoNMeanXnParticle_[c]/(cellVolume*nAvTimeSteps);
                rhoMMean_[c] =  rhoMMeanXnParticle_[c]/(cellVolume*nAvTimeSteps);
                UMean_[c] = momentumMeanXnParticle_[c]/(rhoMMean_[c]*cellVolume*nAvTimeSteps);
                linearKEMean_[c] = 0.5*linearKEMeanXnParticle_[c] / (cellVolume*nAvTimeSteps);                        
                translationalT_[c] = 2.0/(3.0*physicoChemical::k.value()*rhoNMean_[c])*(linearKEMean_[c] - 0.5*rhoMMean_[c]*(UMean_[c] & UMean_[c]));                    
                p_[c] = rhoNMean_[c]*physicoChemical::k.value()*translationalT_[c]; 
            }
            else
            {
                p_[c] = 0.0;
            }   
        }
    }
    
    //*****Ion aceleration calculations***********

    if(time_.samplingTime())
    {
        //Electric field calculation
        scalar cathode1Length = mag(cathode1xi_.x()-cathode1xo_.x());
        scalar cathode2Length = mag(cathode2xi_.x()-cathode2xo_.x());
        //Electric field angles
        scalar theta1; 
        scalar theta2;
        //Distance between the centres of the electric fields in each electrode
        scalar r1;
        scalar r2;
        //Angle between the electrode plates
        scalar E_angle = 90*M_PI/180; 

        //Electron multiplication
        // Correction factor (dimensionless)
        scalar yi; 
        // Correction factor (dimensionless)
        scalar vy;
        // Electron temperature (K)
        scalar Te; 
        // Electron thermal velocity (m/s)
        scalar ue; 
                
        //Constants
        // Electron mass (kg)
        scalar me = 9.11e-31;  
        // Electron charge (C)
        scalar q = 1.602e-19; 
        //Boltzmann constant (J/K)
        scalar K = 1.3806e-23; 
        //Permitivitty in vacuum (F/m)
        scalar epsilon = 8.854e-12;  

        //Variables 
        // Electrode plate length (m) 
        scalar l = cathode1Length;   
        // Electrode plate width (m) 
        scalar w = mag(cathode1xi_.z()-cathode1xo_.z());
        //Electric field
        
        scalar E = voltage_/gap_;
        scalar Q = E*4*M_PI*epsilon*pow(gap_,2.0);
        
        scalar y = 3.79e-5*sqrt(enhancement_*E)/workFunction_;
        scalar v = 1.0 - pow(y,2.0)*(1.0 - 0.333*log(y));
        scalar tSquared = 1.0 + pow(y,2.0)*(1.0 - 0.111*log(y));
        
        //Electron current density
        scalar J = 1.541e-6*pow(enhancement_*E,2.0)/(workFunction_*tSquared)
                    *exp((6.8309*pow(workFunction_,1.5))/(enhancement_*E));
                    
        Info << "J = " << J << endl;
        Info << "I = J*A = " << J*3.025e-15 << endl;
        
        //Electron drift velocity
        scalar u = J/(1e29*q);
        
        Info << "u = " << u << endl;
        
        // Electron density (m^-3)
        //scalar ne = (E*4.0*M_PI*epsilon*pow(gap_,2))/(q*l*w*gap_);
        scalar ne = J/(q*u);

        ionCounter_ = 0.0;

        //cell coordinates range   
        scalar regionx_ = mag(cellzoneStart_.x() - cellzoneEnd_.x());
        scalar regiony_ = mag(cellzoneStart_.y() - cellzoneEnd_.y());

        scalar cellsizex_ = regionx_/x_;   
        scalar cellsizey_ = regiony_/y_;  

        scalar deltax_ = cellsizex_*0.5;
        scalar deltay_ = cellsizey_*0.5;

        scalarField coordinateX_(x_, scalar(0.0)); 
        scalarField coordinateY_(y_, scalar(0.0));

        scalar Eix = 0.0;
        scalar Eiy = 0.0;
        scalar theta = 0.0;
        
        for(label i = 0; i < x_; i++)
        {
            coordinateX_[i] =  cellzoneStart_.x()+ deltax_ + i*cellsizex_;
        }
        
        for(label j = 0; j < y_; j++)
        {
            coordinateY_[j] =  cellzoneStart_.y()+ deltay_ + j*cellsizey_;
        }

        for(label i= 0; i <= x_ - 1; i++)
        {
            // Limit the calculations just for the electrodes' area
            if ((anode1xo_.x() <= coordinateX_[i] && coordinateX_[i] < anode1xi_.x()) || (anode2xo_.x() <= coordinateX_[i] &&  coordinateX_[i] < anode2xi_.x())) 
            {
                //*****Electric field calculations***
                    
                // Electric field in the anode region where there is no parallel cathode
                if(nanode_ == 1 && ncathode_ > 1 && ( cathode1xi_.x() <= coordinateX_[i] && coordinateX_[i] <= cathode2xo_.x()))
                {
                    theta1 = atan(gap_/(0.5*cathode1Length+(coordinateX_[i]-cathode1xi_.x())));
                    theta2 = atan(gap_/(0.5*cathode2Length+(cathode2xo_.x()-coordinateX_[i])));
                    r1 = (0.5*cathode1Length+(coordinateX_[i]-cathode1xi_.x()))/cos(theta1);
                    r2 = (0.5*cathode2Length+(cathode2xo_.x()-coordinateX_[i]))/cos(theta2);
                                
                    Eix = -(Q/(4*M_PI*epsilon*pow(r1,2)))*cos(theta1) + (Q/(4*M_PI*epsilon*pow(r2,2)))*cos(theta2);
                    Eiy = (Q/(4*M_PI*epsilon*pow(r1,2)))*sin(theta1) + (Q/(4*M_PI*epsilon*pow(r2,2)))*sin(theta2);
                    E = sqrt(pow(Eix,2)+pow(Eiy,2));
                }
                
                // Electric field-left side in discontineous parallel plates configuration
                if(ncathode_ > 1  && ((cathode1xo_.x() < coordinateX_[i] && coordinateX_[i] < cathode1xi_.x()) && (anode1xo_.x() < coordinateX_[i] && coordinateX_[i] < anode1xi_.x())))
                {
                    theta2 = atan(gap_/(0.5*cathode2Length+cathode2xo_.x()-coordinateX_[i]));
                    r2=(0.5*cathode2Length+(cathode2xo_.x()-coordinateX_[i]))/cos(theta2);
                    Eix = (voltage_/gap_)*cos(E_angle) + (Q/(4*M_PI*epsilon*pow(r2,2)))*cos(theta2);
                    Eiy = (voltage_/gap_)*sin(E_angle) + (Q/(4*M_PI*epsilon*pow(r2,2)))*sin(theta2);
                    E = sqrt(pow(Eix,2)+pow(Eiy,2));						
                }
                
                // Electric field-right side in discontineous parallel plates configuration        
                if(ncathode_ > 1  && ( cathode2xo_.x() < coordinateX_[i] && coordinateX_[i] < cathode2xi_.x()) && (( anode1xo_.x() < coordinateX_[i] && coordinateX_[i] < anode1xi_.x()) || ( anode2xo_.x() < coordinateX_[i] && coordinateX_[i] < anode2xi_.x())))
                {
                    theta1 = atan(gap_/(0.5*cathode1Length+(cathode2xo_.x()-cathode1xi_.x())+(coordinateX_[i]-cathode2xo_.x())));
                    r1 = (0.5*cathode1Length+(cathode2xo_.x()-cathode1xi_.x())+(coordinateX_[i]-cathode2xo_.x()))/cos(theta1);
                    Eix = (voltage_/gap_)*cos(E_angle) - (Q/(4*M_PI*epsilon*pow(r1,2)))*cos(theta1);
                    Eiy = (voltage_/gap_)*sin(E_angle) + (Q/(4*M_PI*epsilon*pow(r1,2)))*sin(theta1);
                    E = sqrt(pow(Eix,2)+pow(Eiy,2));	
                }
                      
                //Electric field in uniform parallel electrode plates configuration 
                if(nanode_ == 1 && ncathode_  == 1)
                {
                    Eix = (voltage_/gap_)*cos(E_angle);
                    Eiy = (voltage_/gap_)*sin(E_angle);
                    E = sqrt(pow(Eix,2)+pow(Eiy,2));  
                }
               
                //********Electron multiplication*************    
                scalar pressuren = 0.0;
                scalar pressuresum = 0.0;
                scalar pressureavrg = 0.0;
                    
                //--Pressure average calculation for cell columns---- 
                for(label j = 0; j <= y_-1; j++)
                {
                    forAll(controlZone(), c)
                    {
                        const label& cellI = controlZone()[c];
                        const point& cC = cloud_.mesh().cellCentres()[cellI];
                        
                        if(((coordinateX_[i]-deltax_) <= cC.x() &&  cC.x()< (coordinateX_[i]+deltax_)) && ((coordinateY_[j]-deltay_) <= cC.y() && cC.y() < (coordinateY_[j]+deltay_)))
                        {
                            //don't consider the cell in the edges. The concentration of particles by the electric field at edges will modify drastically the pressure average.
                                                            
                            if(j > 1 && j < y_-2) 
                            {
                                pressuren += 1.0;
                                pressuresum += p_[cellI];
                            }
                        } 
                    }
                }
                
                if(Pstream::parRun())
                {
                    reduce(pressuresum, sumOp<scalar>());
                    reduce(pressuren, sumOp<scalar>());
                }
        
                if(pressuren != 0)  
                {
                    pressureavrg = pressuresum/pressuren;
                }
                
                // A constant in function of pressure by fitting the curve using asymmetrical sigmoidal  with experimental data (publications)
                scalar A = 0.0;
                
                // B constant in fucntion of pressure by fitting the curve using asymmetrical sigmoidal
                scalar B = 0.0;
                
                //ionization coefficient (cm^-1)
                scalar alpha = 0.0;
                
                scalar ealphad = 0.0;
                
                // Difference of electron density from cathode to anode (dimensionless)
                scalar dN = 0.0;
            
                //---N2 gas constant for breakdown voltage equation-------------
                if(pressureavrg > VSMALL )
                {   
                    
                    A = 0.83 + (15.5 - 0.83)/(1.0 + pow(pressureavrg/4516,2)); 
                    B = 317.0 + (27735.0 - 317.0)/(1.0 + pow(pressureavrg/1645,2.24));   
                }
                    
                alpha = A*pressureavrg*0.0075*exp(-B*pressureavrg*0.0075/(E/100));   
                ealphad = exp(alpha*gap_*100);
                    
                //------For microgaps > 100 microns------------------
                        
                if(gap_ > 100e-6)
                {
                    dN = ealphad/(1-sEE_*(ealphad-1));
                }
                    
                //----For microgaps < 100 um ----------
        
                if(gap_ <= 100e-6)
                {
                    //correction factor for electron field emission  (dimensionless)
                    yi  = 0.0003795*sqrt(E)/workFunction_;
                            
                    //Message in case the voltage is higher than feasible results
                    if(yi > 1.5)
                    {
                            FatalErrorIn("microDCdischargeAccelerationController::controlParcelsBeforeMove()")
                            << "Provide a lower voltage -> the electric field is stronger to allow calcultions  " << nl << "in: "
                            << time_.time().system()/"controllersDict"
                            << exit(FatalError);
                    }
                        
                    // v(y) correction factor for electron field emission (dimensionless)
                    vy = -0.7519*pow(yi,2)-0.2744*yi + 1.0161;
                
                    //secondary electron emission coefficient (dimensionless)
                    sEEm_ = 1e7*exp(-6.8e7*pow(workFunction_,1.5)*vy/((E/100)*enhancement_)); 
                    
                    //difference of current density from cathode to anode(dimensionless)
                    dN = ealphad/(1-(sEEm_+sEE_)*(ealphad-1));                 
                }

                //-----General calculations-----------------  

                // Electron temperature (K)
                Te = (2.0*epsilon*pow(E,2))/(ne*dN*K); 
                
                Info << "ne = " << ne << endl;
                Info << "dN = " << dN << endl;
                Info << "K = " << K << endl;
                Info << "Te = " << Te << endl;
                
                // Electron velocity average (m/s)
                ue = sqrt(8.0*Te*K/((M_PI)*me)); 
                
                //---------------Maxwell-Boltzmann distribution---------        
                            
//                 scalar EEDlengthControl = 0.0;
                   
//                 //velocity average probability calculation for normalization
//                 scalar EEDaverage =(pow(me/(K*Te),1.5))*4*M_PI*pow(ue,2)*exp(-me*pow(ue,2)/(2*K*Te)); 
//                 
//                 scalar EDDxMax = 0.9999;
//                 
//                 //maximum velocity to set the maximum point in y-axis for normailization of EED (m/s)
//                 scalar uEEDmax =sqrt((-2*K*Te*log(1-EDDxMax))/me); 
// 
//                 // x-axis division of electron velocity for the Maxwell-Boltzmann distribution  function (dimensionless) 
//                 scalar ueEDDxDivision = uEEDmax/(y_ + 2.0);                
//             
//                 //counter for while
//                 label j = 0;
//                 
//                 scalar uDivision = ueEDDxDivision;
//                     
//                 while(EEDlengthControl <= EDDxMax)
//                 {
//                     // Maxwell-Boltzmann distribution velocity function
//                     scalar EED =pow(me/(K*Te),1.5)*4*M_PI*pow(uDivision,2)*exp(-me*pow(uDivision,2)/(2*K*Te));
//                     
//                     scalar EDDNorm = EED/EEDaverage;
//                 
//                     // Function to control the size distribution of the x-axis
//                     EEDlengthControl = 1.0 - exp(-me*pow(uDivision,2)/(2*K*Te));                   
//                                 
//                     //-----Ionic calculations--------- 
//                     if(gap_ <= 15e-6)
//                     {
//                         // Ion density (m^-3) for gaps below 15 microns
//                         rhoI_ = ne*(dN - 1.0)*EDDNorm; 
//                     }
//                     
//                     if(15e-6 < gap_ && gap_ <= 100e-6 && dN >= 1.67) 
//                     {
//                         // Ion density (m^-3) for gaps between 15 and 100 microns where the voltage is enough to reach breakdown voltage
//                         rhoI_ = ne*(dN - 1.67)*EDDNorm;
//                     }
//                     if(15e-6 < gap_ && gap_ <= 100e-6  && dN < 1.67) 
//                     {
//                         // Ion density (m^-3) for gaps between 15 and 100  microns where the voltage is not enough to reach breakdown voltage
//                         rhoI_ = 0;
//                     }
//                             
//                     if(rhoI_> VSMALL) 
//                     {
//                         theta = atan(Eiy/Eix);
//                         
//                         // Set the direction (+/-) of the electric field
//                         if(Eix < VSMALL)
//                         {
//                             theta = (180*M_PI/180)+theta;
//                         }
//                         
//                         forAll(cells, c)
//                         {                
//                             const label& cellI = cells[c];
//                             const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
//                             const point& cC = cloud_.mesh().cellCentres()[cellI];
//                             
//                             current_rhoI_[cellI] = 0;
//                             rhoNMeanXnIonParticle_[cellI] = 0; 
//                             
//                             //normalize temperature constant calculations
//                              //atmospheric pressure reference
//                             scalar pressureRef = 100000;
//                             
//                             scalar TiNormRef =  6800.0 + (29700.0 - 6800.0)/pow(1 + pow(pressureRef/1.5,5.7),10.0);
//                             
//                             scalar TiNorm =  (6800.0 + (29700.0 - 6800.0)/pow(1 + pow(pressureavrg_/1.5,5.7),10.0))/TiNormRef;
//                              
//                             // Temperature of the ion particle (K)
//                             scalar Ti = (epsilon*pow(E,2))/((rhoI_+ne)*K)*TiNorm;
//                                   
//                             dsmcParcel* p = molsInCell[0];
//                             
//                             const dsmcParcel::constantProperties& constProp 
//                                 = cloud_.constProps(p->typeId());
//                             
//                             const scalar& mass = constProp.mass();
//                                                         
//                             // Root mean square ion thermal velocity (m/s)
//                             scalar uIT = sqrt(Ti*K/(M_PI*mass));    
//                             
//                             accelerations_[cellI] = vector::zero;
//                             particlesToAccelerate_[cellI].clear();
//                                         
//                             //-------Particles accelerator-------------
//                             if(((coordinateX_[i]-deltax_) <= cC.x() &&  cC.x() < (coordinateX_[i]+deltax_)) && ((coordinateY_[j]-deltay_) <= cC.y() && cC.y() < (coordinateY_[j]+deltay_)))
//                             {                               
//                                 if(rhoNMean_[cellI] > VSMALL)
//                                 {
//                                     scalar currentRhoI = 0.0;
//                                     scalar ax = 0.0;
//                                     scalar ay = 0.0;
//                                     scalar az = 0.0;
//                                     
//                                     if(Eix < 0)
//                                     {
//                                         // Acceleration for x direction if the electric field vector in x is negative
//                                         ax = -1*pow(uIT*cos(theta) ,2)/gap_; 
//                                     }
//                                                     
//                                     if(Eix >= 0)
//                                     {
//                                         // Acceleration for x direction if the electric field vector in x is positive
//                                         ax = pow(uIT*cos(theta) ,2)/gap_;
//                                     }
//                                     
//                                     // Acceleration for y direction 
//                                     ay = pow(uIT*sin(theta) ,2)/gap_; 
//                                     
//                                     // Acceleration for z direction( is a 2D simulation)
//                                     az = 0; 
//                                     
//                                     accelerations_[cellI].x() = ax;
//                                     accelerations_[cellI].y() = ay;
//                                     accelerations_[cellI].z() = az;
//                                     
//                                     const scalar& nParticle = cloud_.nParticle(); 
//                                     const scalar& cellVolume = mesh_.cellVolumes()[cellI];
//                                     
// //                                     Info << "rhoI_ = " << rhoI_ << endl;
//                                     
//                                     //do
//                                     //{
//                                         forAll(cellOccupancy[cellI], mIC)
//                                         {
//                                             dsmcParcel* p = molsInCell[mIC];
//                                             
//                                             // Control the amount of accelerated particles in each cell
//                                             particlesToAccelerate_[cellI].append(p);
//                                             
//                                             ionCounter_++; 
//                                             
//                                             //count of the total number of accelerated particles(control for accelerate the specific number of ions in each cell)
//                                             rhoNMeanXnIonParticle_[cellI] += nParticle; 
//                                             scalar rhoNMeanXnIonParticle = rhoNMeanXnIonParticle_[cellI];
//                                             currentRhoI =  rhoNMeanXnIonParticle/cellVolume; 
// //                                             Info << "currentRhoI = " << currentRhoI << endl;
//                                             if(currentRhoI >= rhoI_)
//                                             {
//                                                 break;
//                                             }
//                                         }  
//                                     //}
//                                     //while(currentRhoI < rhoI_); 
//                                     
//                                     //Info << "currentRhoI = " << currentRhoI << endl;
//                                     Info << "rhoI_ = " << rhoI_ << endl;
// 
//                                     forAll(particlesToAccelerate_[cellI], mIC)
//                                     {
//                                         dsmcParcel* p = particlesToAccelerate_[cellI][mIC];
//                                         p->U() += 0.5*accelerations_[cellI]*mesh_.time().deltaTValue();
//                                     }          
//                                 }
//                             }
//                         }
//                     }
//                     
//                     // Count the distribution area below the curve that it is being cover until reach the numeric value equivalent to the 99.99%
//                     uDivision += ueEDDxDivision; 
//                       
//                     // Change of cell in y-direction
//                     j++;
//                 }
            }      
        }
    }
    
    scalar ionCounterReduce = ionCounter_;

    if(Pstream::parRun())
    {   
        reduce(ionCounterReduce, sumOp<scalar>());   
    }

    const scalar& nParticle = cloud_.nParticle();
    scalar ionDensity = ionCounterReduce*nParticle; 

    Info <<  " Number of DSMC ion particles    = "  << ionCounterReduce <<nl << endl;
    Info <<  " Total ion density (m^-3)        = "  << ionDensity <<nl << endl;

}
  
void microDCdischargeAccelerationController::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = time_.time();
    if(runTime.outputTime())
    {
    }
}

void microDCdischargeAccelerationController::controlParcelsBeforeCollisions()
{

}

void microDCdischargeAccelerationController::controlParcelsAfterCollisions()
{
    const labelList& cells = mesh_.cellZones()[regionId_];
    
    forAll(cells, c)
    {
        const label& cellI = cells[c];
        
        forAll(particlesToAccelerate_[cellI], mIC)
        {
            dsmcParcel* p = particlesToAccelerate_[cellI][mIC];
            p->U() += 0.5*accelerations_[cellI]*mesh_.time().deltaTValue();
        }
    }
}

void microDCdischargeAccelerationController::updateProperties(const dictionary& newDict)
{
    updateStateControllerProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    setProperties();
}

void microDCdischargeAccelerationController::setProperties()
{
    cellzoneStart_ = propsDict_.lookup("startPoint");
    cellzoneEnd_ = propsDict_.lookup("endPoint");
    x_ = readScalar(propsDict_.lookup("number_of_cells_in_xZone"));   
    y_ = readScalar(propsDict_.lookup("number_of_cells_in_yZone")); 

    //Electric Source properties
    voltage_ = readScalar(propsDict_.lookup("voltage"));                       

    //Micro-gaps breakdown voltage properties 
    workFunction_ = readScalar(propsDict_.lookup("workFunction"));   
    enhancement_ = readScalar(propsDict_.lookup("enhancementFactor")); 
    sEE_ = readScalar(propsDict_.lookup("SEECoeff"));

    //Electrode's geometry
    gap_ = readScalar(propsDict_.lookup("gap"));
    ncathode_ = readScalar(propsDict_.lookup("number_of_cathodes"));
    nanode_ = readScalar(propsDict_.lookup("number_of_anodes"));
    cathode1xo_ = propsDict_.lookup("cathode1_start_point");
    cathode1xi_ = propsDict_.lookup("cathode1_end_point");
    cathode2xo_ = propsDict_.lookup("cathode2_start_point");
    cathode2xi_ = propsDict_.lookup("cathode2_end_point");
    anode1xo_ = propsDict_.lookup("anode1_start_point");
    anode1xi_ = propsDict_.lookup("anode1_end_point");
    anode2xo_ = propsDict_.lookup("anode2_start_point");
    anode2xi_ = propsDict_.lookup("anode2_end_point");   
}
}
// End namespace Foam
// ************************************************************************* //
