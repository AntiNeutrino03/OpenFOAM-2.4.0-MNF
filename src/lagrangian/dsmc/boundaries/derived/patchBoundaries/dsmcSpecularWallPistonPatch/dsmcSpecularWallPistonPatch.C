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

#include "dsmcSpecularWallPistonPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcSpecularWallPistonPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcSpecularWallPistonPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcSpecularWallPistonPatch::dsmcSpecularWallPistonPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
    
    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcSpecularWallPistonPatch::~dsmcSpecularWallPistonPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcSpecularWallPistonPatch::initialConfiguration()
{}

void dsmcSpecularWallPistonPatch::calculateProperties()
{}

void dsmcSpecularWallPistonPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{

    measurePropertiesBeforeControl(p);

    vector& U = p.U();
    
    scalar& pPos = p.position().x();
    
    const scalar& wallPos = p.position().x();
    
    const scalar& deltaT = mesh_.time().deltaTValue();
    
//     scalar trackPercentage = (pPos - p.initialPosition().x())/(pPos + (p.initialPosition().x()*U.x()*deltaT))
//     
    scalar& stepFraction = p.stepFraction();
    
//     Info << "pPos = " << pPos << endl;
//     Info << "p.initialPosition().x() = " << p.initialPosition().x() << endl;

    //scalar pEndPos = pPos + p.initPos().x()*deltaT*U.x();
    
     //scalar initialPos = pPos - stepFraction*deltaT*U.x();
     
    scalar endPos = p.initPos().x() + deltaT*U.x();
    
    scalar newPos = 0.0;
    
    if(endPos > wallPos)
    {
        U.x() = 2.0*velocity_.x() - U.x();
        
        //scalar newPos = 2.0*wallPos - p.initialPosition().x() + U.x()*deltaT;
        newPos = 2.0*wallPos - p.initPos().x() + U.x()*deltaT;
    }
    
    if(newPos > wallPos)
    {
        //particle moves beyond piston, delete
        td.keepParticle = false;
    }
    else
    {
         stepFraction = 1.0;
         pPos = newPos;
    }
    
    measurePropertiesAfterControl(p, 0.0);
}

void dsmcSpecularWallPistonPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void dsmcSpecularWallPistonPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
}

void dsmcSpecularWallPistonPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}



} // End namespace Foam

// ************************************************************************* //
