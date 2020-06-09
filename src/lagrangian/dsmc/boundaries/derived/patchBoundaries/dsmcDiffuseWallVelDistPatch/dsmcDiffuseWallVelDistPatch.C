/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "dsmcDiffuseWallVelDistPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseWallVelDistPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcDiffuseWallVelDistPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseWallVelDistPatch::dsmcDiffuseWallVelDistPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    timeDict_(dict.subDict("timeProperties")),
    time_(t, timeDict_),
    UMean_(vector::zero),
    Ucollected_(vector::zero),
    nParcels_(0),
    binWidth_(readScalar(propsDict_.lookup("binWidth"))), 
    distrX_(binWidth_),
    distrY_(binWidth_),
    distrZ_(binWidth_)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseWallVelDistPatch::~dsmcDiffuseWallVelDistPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDiffuseWallVelDistPatch::updateTimeProperties
(
    const dictionary& newDict
)
{
    timeDict_ = newDict.subDict("timeProperties");

    if (timeDict_.found("resetAtOutput"))
    {
        time_.resetFieldsAtOutput() = Switch(timeDict_.lookup("resetAtOutput"));
    }
}

void dsmcDiffuseWallVelDistPatch::initialConfiguration()
{
    
}

void dsmcDiffuseWallVelDistPatch::calculateProperties()
{

}

void dsmcDiffuseWallVelDistPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    const Time& runTime = time_.time();
    
    measurePropertiesBeforeControl(p);
    
//     scalar currentTime = cloud_.mesh().time().value();

    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    labelList& vibLevel = p.vibLevel();
    
    label& ELevel = p.ELevel();

    label typeId = p.typeId();

    vector nw = p.faceNormal();
    nw /= mag(nw);
    
//     if(nw.y() == 1 || nw.y() == -1)
//     {
        //Particle hitting lateral face
        nParcels_++;
        Ucollected_ += p.U(); 
        
        vector Ucollected = Ucollected_;
        label nParcels = nParcels_;

//         //- sending
//         if(Pstream::parRun())
//         {
//             for (int p = 0; p < Pstream::nProcs(); p++)
//             {
//                 if(p != Pstream::myProcNo())
//                 {
//                     const int proc = p;
//                     {
//                         OPstream toNeighbour(Pstream::blocking, proc);
//                         toNeighbour << Ucollected << nParcels;
//                     }
//                 }
//             }
//             
//             Pout << "1" << endl;
//         
//             //- receiving
//             for (int p = 0; p < Pstream::nProcs(); p++)
//             {
//                 if(p != Pstream::myProcNo())
//                 {
//                     vector UcollectedProc;
//                     label nParcelsProc;
//     
//                     const int proc = p;
//                     {
//                         IPstream fromNeighbour(Pstream::blocking, proc);
//                         fromNeighbour >> UcollectedProc >> nParcelsProc;
//                     }
//     
//                     Ucollected += UcollectedProc;
//                     nParcels += nParcelsProc;
//                 }
//             }
//         }
        
        if(nParcels_ > 0)
        {
            UMean_ = Ucollected/nParcels;
        }
        
        distrX_.add((p.U().x() - UMean_.x()));
        distrY_.add((p.U().y() - UMean_.y()));
        distrZ_.add((p.U().z() - UMean_.z()));

        if(runTime.outputTime())
        {
            fileName timePath(runTime.path()/runTime.timeName()/"uniform");
        
            if (!isDir(timePath))
            {
                mkDir(timePath);
            }

            List< Pair<scalar> > rawDistributionX = distrX_.raw();
            List< Pair<scalar> > rawDistributionY = distrY_.raw();
            List< Pair<scalar> > rawDistributionZ = distrZ_.raw();

            label nSizeX = rawDistributionX.size();
            label nSizeY = rawDistributionY.size();
            label nSizeZ = rawDistributionZ.size();
            
//             Pout << "3" << endl;

//             if (Pstream::parRun())
//             {
//                 reduce(nSizeX, sumOp<label>());
//                 reduce(nSizeY, sumOp<label>());
//                 reduce(nSizeZ, sumOp<label>());
//             }
            
//             Pout << "4" << endl;

            scalarField xAxisX (nSizeX, 0.0);
            scalarField yAxisX (nSizeX, 0.0);
            scalarField xAxisY (nSizeY, 0.0);
            scalarField yAxisY (nSizeY, 0.0);
            scalarField xAxisZ (nSizeZ, 0.0);
            scalarField yAxisZ (nSizeZ, 0.0);
            

            forAll(rawDistributionX, i)
            {
                xAxisX[i] = rawDistributionX[i].first();
                yAxisX[i] = rawDistributionX[i].second();
            }
            
            forAll(rawDistributionY, i)
            {
                xAxisY[i] = rawDistributionY[i].first();
                yAxisY[i] = rawDistributionY[i].second();
            }
            
            forAll(rawDistributionZ, i)
            {
                xAxisZ[i] = rawDistributionZ[i].first();
                yAxisZ[i] = rawDistributionZ[i].second();
            }

//             //- sending
//             if (Pstream::parRun())
//             {
//                 for (int p = 0; p < Pstream::nProcs(); p++)
//                 {
//                     if(p != Pstream::myProcNo())
//                     {
//                         const int proc = p;
//                         {
//                             OPstream toNeighbour(Pstream::blocking, proc);
//                             toNeighbour << xAxisX << yAxisX << xAxisY 
//                                         << yAxisY << xAxisZ << yAxisZ;
//                         }
//                     }
//                 }
//             
//                 //- receiving
//                 for (int p = 0; p < Pstream::nProcs(); p++)
//                 {
//                     if(p != Pstream::myProcNo())
//                     {
//                         scalarField xAxisXProc;
//                         scalarField yAxisXProc;
//                         scalarField xAxisYProc;
//                         scalarField yAxisYProc;
//                         scalarField xAxisZProc;
//                         scalarField yAxisZProc;
//         
//                         const int proc = p;
//                         {
//                             IPstream fromNeighbour(Pstream::blocking, proc);
//                             fromNeighbour >> xAxisXProc >> yAxisXProc >> xAxisYProc 
//                                         >> yAxisYProc >> xAxisZProc >> yAxisZProc;
//                         }
//         
//                         forAll(xAxisXProc, i)
//                         {
//                             xAxisX[i] += xAxisXProc[i];
//                             yAxisX[i] += yAxisXProc[i];
//                         }
//                         
//                         forAll(xAxisYProc, i)
//                         {
//                             xAxisY[i] += xAxisYProc[i];
//                             yAxisY[i] += yAxisYProc[i];
//                         }
//                         
//                         forAll(xAxisZProc, i)
//                         {
//                             xAxisZ[i] += xAxisZProc[i];
//                             yAxisZ[i] += yAxisZProc[i];
//                         }
//                     }
//                 }
//             }
            
            word binsX = "velocityBinsX";
            word binsY = "velocityBinsY";
            word binsZ = "velocityBinsZ";
            word X = "velocityDistributionX";
            word Y = "velocityDistributionY";
            word Z = "velocityDistributionZ";

            writeTimeData(timePath, binsX, xAxisX);
            writeTimeData(timePath, binsY, xAxisY);
            writeTimeData(timePath, binsZ, xAxisZ);
            writeTimeData(timePath, X, yAxisX);
            writeTimeData(timePath, Y, yAxisY);
            writeTimeData(timePath, Z, yAxisZ);
        }
//     }

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    Random& rndGen(cloud_.rndGen());

    while (mag(Ut) < SMALL)
    {
        // If the incident velocity is parallel to the face normal, no
        // tangential direction can be chosen.  Add a perturbation to the
        // incoming velocity and recalculate.

        U = vector
        (
            U.x()*(0.8 + 0.2*rndGen.scalar01()),
            U.y()*(0.8 + 0.2*rndGen.scalar01()),
            U.z()*(0.8 + 0.2*rndGen.scalar01())
        );

        U_dot_nw = U & nw;

        Ut = U - U_dot_nw*nw;
    }

    // Wall tangential unit vector
    vector tw1 = Ut/mag(Ut);

    // Other tangential unit vector
    vector tw2 = nw^tw1;

    const scalar& T = temperature_;

    scalar mass = cloud_.constProps(typeId).mass();

    scalar rotationalDof = cloud_.constProps(typeId).rotationalDegreesOfFreedom();
    
    scalar vibrationalDof = cloud_.constProps(typeId).vibrationalDegreesOfFreedom();

    U =
        sqrt(physicoChemical::k.value()*T/mass)
       *(
            rndGen.GaussNormal()*tw1
          + rndGen.GaussNormal()*tw2
          - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
        );
       
    ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);
    
    vibLevel = cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);
    
    ELevel = cloud_.equipartitionElectronicLevel
                    (
                        T,
                        cloud_.constProps(typeId).degeneracyList(),
                        cloud_.constProps(typeId).electronicEnergyList(),
                        typeId
                    );   
    
    U += velocity_;
    
    measurePropertiesAfterControl(p, 0.0);
}

void dsmcDiffuseWallVelDistPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}


void dsmcDiffuseWallVelDistPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}

void dsmcDiffuseWallVelDistPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
}

} // End namespace Foam

// ************************************************************************* //
