Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "createPhi.H"

singlePhaseTransportModel laminarTransport(U, phi);

// Solutal expansion coefficient [m3/mol]
dimensionedScalar betaC(laminarTransport.lookup("betaC"));

// Reference temperature [K]
dimensionedScalar TRef(laminarTransport.lookup("TRef"));

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Set the turbulence into the reaction model
//reaction->setTurbulence(turbulence());

//Info<< "Creating field dpdt\n" << endl;
//volScalarField dpdt
//(
//    IOobject
//    (
//        "dpdt",
//        runTime.timeName(),
//        mesh
//    ),
//    mesh,
//    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
//);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

// Create chemistry model
autoPtr<concChemistryModel> pChemistry // OF 2.3
(
    concChemistryModel::New(mesh) // OF 2.3
);
concChemistryModel& chemistry = pChemistry(); // OF 2.3

// Create thermo model
concReactionThermo& thermo = chemistry.thermo(); // OF 2.3
thermo.validate(args.executable(), "h", "e");

concBasicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();
PtrList<volScalarField>& C = composition.C();

word inertSpecie(thermo.lookup("inertSpecie"));

//volScalarField& p = thermo.p();

volScalarField& hs = thermo.he();  

volScalarField& T = thermo.T();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, speciesI)
{
    fields.add(Y[speciesI]);
}
fields.add(hs);

DimensionedField<scalar, volMesh> chemistrySh
(
    IOobject
    (
        "chemistry::Sh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
);

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);

volScalarField gamma 
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// Reference physical constants
const dimensionedScalar& R = constant::physicoChemical::R;
const dimensionedScalar& F = constant::physicoChemical::F;
//const dimensionedScalar& kB = constant::physicoChemical::k;
//const dimensionedScalar& e = constant::electromagnetic::e;


// Create electrochemistry model
electrochemistryModels::electrochemistryModel electrochemistry(thermo);
PtrList<electrochemistryModels::electrodeModels::electrodeModel>& 
    electrodes(electrochemistry.electrodes());

electrochemistryModels::electrolyteModels::electrolyteModel& 
    electrolyte(electrochemistry.electrolyte());

volScalarField div_il 
(
    IOobject
    (

        "div_il",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("div_il", dimCurrent/dimVolume, 0.0)
);

volScalarField& phiEl = electrolyte.phiE();

volVectorField& id = electrolyte.iD();
volVectorField& il = electrolyte.i();

volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "E",
        dimensionSet(1,1,-3,0,0,-1,0),
        vector(0.0, 0.0, 0.0)
    )
);

// Ionic conductivity of electrolyte
const volScalarField& kappa(electrolyte.kappa());

// Species charge numbers
const PtrList<label>& z = electrolyte.z();

// Species transference numbers
const PtrList<volScalarField>& t = electrolyte.transNumberModel().t();

// Species molar mass
PtrList<dimensionedScalar> M(C.size());

// Species flux vectors
PtrList<volVectorField>& N = electrolyte.N();

// Species diffusivities
const PtrList<volScalarField>& D(electrolyte.D());

// Effective diffusivities
PtrList<volScalarField> Deff(D.size());
forAll(Y, speciesI)
{
    Deff.set
    (
        speciesI,
        new volScalarField 
        (
            IOobject
            (
                "Deff_"+Y[speciesI].name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            D[speciesI]
        )
    );
}

// Species mobilities
PtrList<volScalarField> mu(D.size());

forAll(Y, speciesI)
{
    mu.set
    (
        speciesI,
        new volScalarField(Deff[speciesI]/(R*T))
    );
}

// Species molar weight in kg/mol
forAll(Y, speciesI)
{
    M.set
    (
        speciesI,
        new dimensionedScalar
        (
           "M_"+Y[speciesI].name(),
           dimensionSet(dimMass/dimMoles),
           thermo.composition().W(speciesI)/1000.0
        )
    );
}

volScalarField Dm 
(
    IOobject
    (
        "Dm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    F*z[0]*Deff[0]
);

volScalarField Di
(
    IOobject
    (
        "Di",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Deff[0]
);
// Number of charged species
label nChargedSpecies = 0;
forAll(Y, speciesI)
{
    if (z[speciesI] != 0)
    {
        nChargedSpecies++;
    }
}

// Charge density
volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoE", dimCurrent*dimTime/dimVol, 0.0)
);

// Sum of charge number-concentration-products of each species excluding the 
// inactive charged species that is used to balance for electroneutrality
volScalarField zC 
(
    IOobject
    (
        "zC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zC", dimMoles/dimVol, 0.0)
);

// Relative permittivity [F/m]
//const dimensionedScalar& epsilonR = electrolyte.epsilon();

// Vacuum permittivity [C]
//const dimensionedScalar& epsilon0 = constant::electromagnetic::epsilon0;

//scalar iAvgTime = 0.0;
//scalar phiEsAvgTime = 0.0;
//label timeStepCounter = 0;
//label startCounter = 200;

const IOdictionary& electrochemDict = electrochemistry.dict();
const word transportSolver = electrochemDict.lookup("transportSolver");

// Reference concentration
const dimensionedScalar& CRef = electrolyte.Cstd();


// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.0 - betaC*(C[0] - CRef)
);


Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", g & mesh.C());
surfaceScalarField ghf("ghf", g & mesh.Cf());

p += rhok*gh;

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    mesh.solutionDict().subDict("PIMPLE"),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}
