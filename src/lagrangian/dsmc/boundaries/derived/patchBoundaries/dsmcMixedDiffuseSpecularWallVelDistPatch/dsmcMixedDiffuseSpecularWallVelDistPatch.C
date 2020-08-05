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

#include "dsmcMixedDiffuseSpecularWallVelDistPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcMixedDiffuseSpecularWallVelDistPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcMixedDiffuseSpecularWallVelDistPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcMixedDiffuseSpecularWallVelDistPatch::dsmcMixedDiffuseSpecularWallVelDistPatch
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
    typeIds_(),
    binWidth_(readScalar(propsDict_.lookup("binWidth"))), 
    diffuseFraction_ (readScalar(propsDict_.lookup("diffuseFraction"))),
    distrX_(),
    distrY_(),
    distrZ_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcMixedDiffuseSpecularWallVelDistPatch::~dsmcMixedDiffuseSpecularWallVelDistPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcMixedDiffuseSpecularWallVelDistPatch::updateTimeProperties
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

void dsmcMixedDiffuseSpecularWallVelDistPatch::initialConfiguration()
{
    
}

void dsmcMixedDiffuseSpecularWallVelDistPatch::calculateProperties()
{

}

void dsmcMixedDiffuseSpecularWallVelDistPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
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
    
    label distSize = typeIds_.size();
    
    //Individual species
    distrX_[typeId].add((p.U().x()));
    distrY_[typeId].add((p.U().y()));
    distrZ_[typeId].add((p.U().z()));
    
    //Mixture
    distrX_[distSize].add((p.U().x()));
    distrY_[distSize].add((p.U().y()));
    distrZ_[distSize].add((p.U().z()));

    if(runTime.outputTime())
    {
        fileName timePath(runTime.path()/runTime.timeName()/"uniform");
    
        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        List< List< Pair<scalar> > > rawDistributionX(typeIds_.size() + 1);
        List< List< Pair<scalar> > > rawDistributionY(typeIds_.size() + 1);
        List< List< Pair<scalar> > > rawDistributionZ(typeIds_.size() + 1);
        
        forAll(rawDistributionX, i)
        {
            rawDistributionX[i] = distrX_[i].raw();
            rawDistributionY[i] = distrY_[i].raw();
            rawDistributionZ[i] = distrZ_[i].raw();
        }

        labelList nSizeX(typeIds_.size() + 1);
        labelList nSizeY(typeIds_.size() + 1);
        labelList nSizeZ(typeIds_.size() + 1);
        
        forAll(nSizeX, i)
        {
            nSizeX[i] = rawDistributionX[i].size();
            nSizeY[i] = rawDistributionY[i].size();
            nSizeZ[i] = rawDistributionZ[i].size();
        }            

        List <scalarField> xAxisX (typeIds_.size() + 1);
        List <scalarField> yAxisX (typeIds_.size() + 1);
        List <scalarField> xAxisY (typeIds_.size() + 1);
        List <scalarField> yAxisY (typeIds_.size() + 1);
        List <scalarField> xAxisZ (typeIds_.size() + 1);
        List <scalarField> yAxisZ (typeIds_.size() + 1);
        
        forAll(xAxisX, i)
        {
            xAxisX[i].setSize(nSizeX[i]);
            yAxisX[i].setSize(nSizeX[i]);
            xAxisY[i].setSize(nSizeY[i]);
            yAxisY[i].setSize(nSizeY[i]);
            xAxisZ[i].setSize(nSizeZ[i]);
            yAxisZ[i].setSize(nSizeZ[i]);
        }
        

        forAll(rawDistributionX, i)
        {
            forAll(rawDistributionX[i], j)
            {
                xAxisX[i][j] = rawDistributionX[i][j].first();
                yAxisX[i][j] = rawDistributionX[i][j].second();
            }
        }
        
        forAll(rawDistributionY, i)
        {
            forAll(rawDistributionY[i], j)
            {
                xAxisY[i][j] = rawDistributionY[i][j].first();
                yAxisY[i][j] = rawDistributionY[i][j].second();
            }
        }
        
        forAll(rawDistributionZ, i)
        {
            forAll(rawDistributionZ[i], j)
            {
                xAxisZ[i][j] = rawDistributionZ[i][j].first();
                yAxisZ[i][j] = rawDistributionZ[i][j].second();
            }
        }

        
        wordList binsX(typeIds_.size() + 1);
        wordList binsY(typeIds_.size() + 1);
        wordList binsZ(typeIds_.size() + 1);
        wordList X(typeIds_.size() + 1);
        wordList Y(typeIds_.size() + 1);
        wordList Z(typeIds_.size() + 1);
        
        const List<word> molecules (propsDict_.lookup("typeIds"));

        DynamicList<word> moleculesReduced(0);

        forAll(molecules, i)
        {
            const word moleculeName(molecules[i]);

            if(findIndex(moleculesReduced, moleculeName) == -1)
            {
                moleculesReduced.append(moleculeName);
            }
        }

        moleculesReduced.shrink();
        
        for(label i = 0; i < typeIds_.size(); i++)
        {
            binsX[i] = "velocityBinsX_"+moleculesReduced[i];
            binsY[i] = "velocityBinsY_"+moleculesReduced[i];
            binsZ[i] = "velocityBinsZ_"+moleculesReduced[i];
            X[i] = "velocityDistributionX_"+moleculesReduced[i];
            Y[i] = "velocityDistributionY_"+moleculesReduced[i];
            Z[i] = "velocityDistributionZ_"+moleculesReduced[i];
        }
        
        binsX[typeIds_.size()] = "velocityBinsX_mixture";
        binsY[typeIds_.size()] = "velocityBinsY_mixture";
        binsZ[typeIds_.size()] = "velocityBinsZ_mixture";
        X[typeIds_.size()] = "velocityDistributionX_mixture";
        Y[typeIds_.size()] = "velocityDistributionY_mixture";
        Z[typeIds_.size()] = "velocityDistributionZ_mixture";

        forAll(binsX, i)
        {
            writeTimeData(timePath, binsX[i], xAxisX[i]);
            writeTimeData(timePath, binsY[i], xAxisY[i]);
            writeTimeData(timePath, binsZ[i], xAxisZ[i]);
            writeTimeData(timePath, X[i], yAxisX[i]);
            writeTimeData(timePath, Y[i], yAxisY[i]);
            writeTimeData(timePath, Z[i], yAxisZ[i]);
        }

    }

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    Random& rndGen(cloud_.rndGen());
    
    if (diffuseFraction_ > rndGen.scalar01())
    {
        // Diffuse reflection

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
    }
    else
    {
        // Specular reflection

        if (U_dot_nw > 0.0)
        {
            U -= 2.0*U_dot_nw*nw;
        }
    }
    
    measurePropertiesAfterControl(p, 0.0);
}

void dsmcMixedDiffuseSpecularWallVelDistPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}


void dsmcMixedDiffuseSpecularWallVelDistPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}

void dsmcMixedDiffuseSpecularWallVelDistPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
    
    //  read in the type ids

    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcFreeStreamInflowPatch::dsmcFreeStreamInflowPatch()")
            << "Cannot have zero typeIds being inserd." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    //  set the type ids

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcDiffuseWallVelDistPatch::dsmcDiffuseWallVelDistPatch()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    
    distrX_.setSize(typeIds_.size() + 1);
    
    forAll(distrX_, i)
    {
        distrX_[i](binWidth_);
    }
    
    distrY_.setSize(typeIds_.size() + 1);
    
    forAll(distrY_, i)
    {
        distrY_[i](binWidth_);
    }
    
    distrZ_.setSize(typeIds_.size() + 1);
    
    forAll(distrZ_, i)
    {
        distrZ_[i](binWidth_);
    }
}

} // End namespace Foam

// ************************************************************************* //
