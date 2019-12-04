// Interpolation für vecWeighted bei j=-sigma*grad(vecWeighted)

#ifndef vecWeighted_H
#define vecWeighted_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "coupledFvPatch.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class vecWeighted Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class vecWeighted
:
    public surfaceInterpolationScheme<Type>
{
    // Private data

        // Reference to fvMesh
        const fvMesh& mesh_;

	const volTensorField& sigma_;

        // Demand-driven data

        mutable surfaceScalarField* oDelta_; // distance face - owner cell
        mutable surfaceScalarField* nDelta_; // distance face - neighbour cell center

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const vecWeighted&);

        //- Construct central-differencing weighting factors
        void makeDeltas() const;

protected:

        //- Clear all geometry and addressing
        void clearOut();

public:

    //- Runtime type information
        TypeName("vecWeighted");


    // Constructors

        //- Construct from Istream
        vecWeighted(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
 	    mesh_(mesh),
            sigma_(this->mesh().objectRegistry::template
                lookupObject<volTensorField>(word(is))),
            oDelta_(NULL),
            nDelta_(NULL)
        {}

        //- Construct from faceFlux and Istream
        vecWeighted
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            mesh_(mesh),
            sigma_(this->mesh().objectRegistry::template
                lookupObject<volTensorField>(word(is))),
            oDelta_(NULL),
            nDelta_(NULL)
        {}

	~vecWeighted()
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
