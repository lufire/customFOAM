{
    fvScalarMatrix hsEqn
    (
          fvm::ddt(voidfraction, hs) - fvm::Sp(fvc::ddt(voidfraction), hs)
        + mvConvection->fvmDiv(phi, hs) - fvm::Sp(fvc::div(phi), hs)
	    + fvc::ddt(voidfraction, K) + fvc::div(phi, K)
        + (
	        hs.name() == "e"
	        ? 
            fvc::div
	        (
		        fvc::absolute(phi/fvc::interpolate(rho), U),
		        p,
	            "div(phiv,p)"
	        )
	        : -dpdt
	      )
        - fvm::laplacian(thermo.alpha()/rho*pow(voidfraction,1.5), hs)
        //+ mstm().multiSpeciesHeatSource()
        ==
          chemistrySh/rho
        + fvOptions(voidfraction, hs)
    );


    hsEqn.relax();

    fvOptions.constrain(hsEqn);

    hsEqn.solve();

    fvOptions.correct(hs);

    thermo.correct();

    // A.Alexiou 2015 - Debugging
    //Info << "Check T" << thermo.T();

    Info<< "T min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
