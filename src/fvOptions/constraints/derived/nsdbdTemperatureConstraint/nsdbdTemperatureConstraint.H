/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

Class
    Foam::fv::nsdbdTemperatureConstraint

Group
    grpFvOptionsConstraints

Description
    Fixed temperature equation constraint to replicate a NS-DBD plasma
    actuator operations

Usage
    \verbatim
    nsdbdTemperature
    {
        type            nsdbdTemperatureConstraint;
        active          true;
        
        nsdbdTemperatureConstraintCoeffs
        {
            selectionMode       cellSet;
            cellSet             set1;
            ambientTemperature  288.15;
            pulseTime           100e-9;
            restTime            0.0099999;//every 0.01 seconds i.e. 100Hz
            actuatorLength      0.006;
            actuatorStartPoint  0.230;
            chi                 0.7;
        }
        nsdbdTemperatureConstraintCoeffs
        {
            selectionMode       cellSet;
            cellSet             set2;
            ambientTemperature  293;
            pulseTime           100e-9;
            restTime            4e-3;
            actuatorLength      0.005;
            actuatorStartPoint  0.025;
            chi                 0.5;
        }
    }
    \endverbatim

Note:
    See D. Gaitonde, "Analysis of plasma-based flow control mechanisms through large-eddy simulations", Computers \& Fluids, 85(1), 19-26, 2013.

SourceFiles
    fvOption.C

\*---------------------------------------------------------------------------*/

#ifndef nsdbdTemperatureConstraint_H
#define nsdbdTemperatureConstraint_H

#include "fvOption.H"
#include "NamedEnum.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class nsdbdTemperatureConstraint Declaration
\*---------------------------------------------------------------------------*/

class nsdbdTemperatureConstraint
:
    public option
{

public:


protected:

    // Protected data

        //- Ambient temperature [K]
        autoPtr<DataEntry<scalar> > Tambient_;
        
        scalar lastOnOffTime_;
        
        //- How long a single actuator pulse lasts
        scalar pulseTime_;
        
        //- Time between pulses
        scalar restTime_;
        
        //- How long the actuator is
        scalar actuatorLength_;
        
        //- x co-ordinate in domain where actuator starts
        scalar actuatorStartPoint_;
        
        //- Parameter for setting the average temperature in the
        // deposition volume (if chi_ = 1.0, the average temperature 
        // is double the ambient temperature)
        scalar chi_;
               
        //- Keep track of whether actuator is resting or pulsing
        bool pulsing_;
        
        bool resting_;
        


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        nsdbdTemperatureConstraint(const nsdbdTemperatureConstraint&);

        //- Disallow default bitwise assignment
        void operator=(const nsdbdTemperatureConstraint&);


public:

    //- Runtime type information
    TypeName("nsdbdTemperatureConstraint");


    // Constructors

        //- Construct from components
        nsdbdTemperatureConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~nsdbdTemperatureConstraint()
    {}


    // Member Functions

        //- Constrain energy equation to fix the temperature
        virtual void constrain(fvMatrix<scalar>& eqn, const label fieldi);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
