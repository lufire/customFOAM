// Interpolation für volWeighted bei j=-sigma*grad(volWeighted)

#ifndef volWeighted_H
#define volWeighted_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class volWeighted Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class volWeighted
:
    public surfaceInterpolationScheme<Type>
{
    // Private data

        // Reference to fvMesh
        const fvMesh& mesh_;

	const volScalarField& sigma_;

        // Demand-driven data

        mutable surfaceScalarField* oDelta_; // distance face - owner cell
        mutable surfaceScalarField* nDelta_; // distance face - neighbour cell center

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const volWeighted&);

        //- Construct central-differencing weighting factors
        void makeDeltas() const;

protected:

        //- Clear all geometry and addressing
        void clearOut();

public:

    //- Runtime type information
        TypeName("volWeighted");


    // Constructors

        //- Construct from Istream
        volWeighted(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
 	    mesh_(mesh),
            sigma_(this->mesh().objectRegistry::template
                lookupObject<volScalarField>(word(is))),
            oDelta_(NULL),
            nDelta_(NULL)
        {}

        //- Construct from faceFlux and Istream
        volWeighted
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            mesh_(mesh),
            sigma_(this->mesh().objectRegistry::template
                lookupObject<volScalarField>(word(is))),
            oDelta_(NULL),
            nDelta_(NULL)
        {}

	~volWeighted()
	{
    	    clearOut();
	}


    // Member Functions

    const surfaceScalarField& oDelta() const;
    const surfaceScalarField& nDelta() const;

    // interpolates cell values to faces
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
    interpolate
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf // volume field
    ) const;


    //- Return the interpolation weighting factors
    tmp<surfaceScalarField> weights
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    ) const
    {
        return this->mesh().surfaceInterpolation::weights();
    }

};
} // End namespace Foam


#endif
