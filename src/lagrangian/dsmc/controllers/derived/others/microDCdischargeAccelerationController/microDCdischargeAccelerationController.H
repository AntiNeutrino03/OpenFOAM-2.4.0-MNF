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

\*---------------------------------------------------------------------------*/

#ifndef microDCdischargeAccelerationController_H
#define microDCdischargeAccelerationController_H

#include "dsmcStateController.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
Class microDCdischargeAccelerationController Declaration
\*---------------------------------------------------------------------------*/

class microDCdischargeAccelerationController
:
        public dsmcStateController
{

private:
        //- subDictionary containing the properties
        dictionary propsDict_;
        labelList typeIds_;
        
        // Pressure measurement
        scalar nTimeSteps_;
        scalarField rhoNMean_;
        scalarField rhoNMeanXnParticle_;
        scalarField rhoMMeanXnParticle_;
        vectorField momentumMeanXnParticle_;
        scalarField linearKEMeanXnParticle_;

        scalarField rhoN_;
        scalarField rhoNMean2_;
        scalarField rhoMMean_;
        scalarField linearKEMean_;
        vectorField UMean_;
        scalarField translationalT_;
        scalarField p_;
        scalar densityOnly_;
        scalar controlVolume_;
        label measuredParcels_;
        scalar avParcelDensity_;
        scalar pressuren_;
        scalar pressuresum_;
        scalar pressureavrg_;
        
        //List of particles to accelerate
        List<DynamicList<dsmcParcel*> > particlesToAccelerate_;
        //Accelerations per cell
        List<vector>  accelerations_;
        
        //Cell ordering by (i,j)
        vector cellzoneStart_;
        vector cellzoneEnd_;
        scalar x_;
        scalar y_;

                            
        //Electric Source properties
        scalar voltage_;

        //work function of the electrode's material (eV)
        scalar workFunction_; 
        //Electric field enhacement factor for micro-gap
        scalar enhancement_;
        //Secondary electron emission coefficient
        scalar sEE_;
        //Secondary electron emission coefficient for microgaps
        scalar sEEm_;
        //Ion number density
        scalar rhoI_;
        scalarField current_rhoI_;
        scalarField  rhoNMeanXnIonParticle_;
        scalarField  rhoNMeanXnIonParticle2_;

        //Electrode's geometry
        scalar gap_;
        scalar ncathode_;
        scalar nanode_;
        vector cathode1xo_;
        vector cathode1xi_;
        vector cathode2xo_;
        vector cathode2xi_;
        vector anode1xo_;
        vector anode1xi_;
        vector anode2xo_;
        vector anode2xi_;

        scalar ionCounter_;
        
        word fieldName_;
                    
        // Private Member Functions
        void setProperties();

    public:

        //- Runtime type information
        TypeName("microDCdischargeAccelerationController");

        // Constructors
        //- Construct from components
        microDCdischargeAccelerationController
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );

        // Destructor
        ~microDCdischargeAccelerationController();

        // Member Functions
        void initialConfiguration();
        void controlParcelsBeforeMove();
        void controlParcelsBeforeCollisions();
        void controlParcelsAfterCollisions();

        //- calculate any required properties
        void calculateProperties();

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
