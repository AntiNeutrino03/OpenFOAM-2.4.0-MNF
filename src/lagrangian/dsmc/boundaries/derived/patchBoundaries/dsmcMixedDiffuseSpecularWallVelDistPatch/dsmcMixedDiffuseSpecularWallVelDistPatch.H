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

Class
    dsmcMixedDiffuseSpecularWallVelDistPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    dsmcMixedDiffuseSpecularWallVelDistPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcMixedDiffuseSpecularWallVelDistPatch_H
#define dsmcMixedDiffuseSpecularWallVelDistPatch_H

#include "dsmcPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"
#include "writeTimeData.H"
#include "timeDataMeas.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcMixedDiffuseSpecularWallVelDistPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcMixedDiffuseSpecularWallVelDistPatch
:
	public dsmcPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;
        
        dictionary timeDict_;
        timeDataMeas time_;
        labelList typeIds_;      

        scalar binWidth_;
         //- Fraction of wall interactions that are diffuse
        scalar diffuseFraction_;
        List <distribution> distrX_;
        List <distribution> distrY_;
        List <distribution> distrZ_;
        
        void setProperties();


public:

    //- Runtime type information
    TypeName("dsmcMixedDiffuseSpecularWallVelDistPatch");

    // Constructors

        //- Construct from components
        dsmcMixedDiffuseSpecularWallVelDistPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcMixedDiffuseSpecularWallVelDistPatch();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        void controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td);

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );
        
        void writeField();

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);
        
        void updateTimeProperties(const dictionary&);



        // Access





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
