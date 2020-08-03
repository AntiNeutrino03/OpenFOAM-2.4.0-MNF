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

\*---------------------------------------------------------------------------*/

#include "nsdbdTemperatureConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(nsdbdTemperatureConstraint, 0);
        addToRunTimeSelectionTable
        (
            option,
            nsdbdTemperatureConstraint,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::nsdbdTemperatureConstraint::nsdbdTemperatureConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    Tambient_(NULL),
    lastOnOffTime_(0.0),
    pulseTime_(readScalar(coeffs_.lookup("pulseTime"))),
    restTime_(readScalar(coeffs_.lookup("restTime"))),
    actuatorLength_(readScalar(coeffs_.lookup("actuatorLength"))),
    actuatorStartPoint_(readScalar(coeffs_.lookup("actuatorStartPoint"))),
    chi_(readScalar(coeffs_.lookup("chi"))),
    pulsing_(true),
    resting_(false)
{      
    Tambient_.reset
    (
        DataEntry<scalar>::New("ambientTemperature", coeffs_).ptr()
    );

    // Set the field name to that of the energy field from which the temperature
    // is obtained

    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.setSize(1, thermo.he().name());

    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::nsdbdTemperatureConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);
    
    scalar currentTime = mesh().time().value();
    
    if(pulsing_ && (currentTime - lastOnOffTime_) <= (pulseTime_ + SMALL))
    {
        Info << "Setting temperature" << endl;
        scalarField Tuni(cells_.size(), Tambient_->value(currentTime));
        
        forAll(cells_, i)
        {
            scalar cC = mesh_.C()[cells_[i]].x();
            
            scalar nonDimPos = (cC-actuatorStartPoint_)/actuatorLength_;
            
            Tuni[i] =  Tambient_->value(currentTime)
                        + chi_*(Tambient_->value(currentTime)
                        *sqrt(1.0/(2.0*3.1416*pow(nonDimPos,3)))
                        *exp(-(1.0/(2.0*0.3*0.3*nonDimPos))
                            *sqr(nonDimPos - 0.3))
                        );
            
            if(Tuni[i] < Tambient_->value(currentTime))
            {
                Tuni[i] = Tambient_->value(currentTime);
            }
        }
        
        eqn.setValues(cells_, thermo.he(thermo.p(), Tuni, cells_));
    }
    
    if(pulsing_ && (currentTime - lastOnOffTime_) > pulseTime_)
    {
        pulsing_ = false;
        resting_ = true;
        lastOnOffTime_ = currentTime;
    }
    
    if(resting_ && (currentTime - lastOnOffTime_) >= (restTime_ + SMALL))
    {
        pulsing_ = true;
        resting_ = false;
        lastOnOffTime_ = currentTime;
    }
}


bool Foam::fv::nsdbdTemperatureConstraint::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        if (coeffs_.found(Tambient_->name()))
        {
            Tambient_.reset
            (
                DataEntry<scalar>::New(Tambient_->name(), dict).ptr()
            );
        }
        
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
