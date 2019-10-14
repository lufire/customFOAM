{
    fvScalarMatrix hsEqn
    (
          fvm::ddt(rho, hs) + mvConvection->fvmDiv(phi, hs)
	    + fvc::ddt(rho, K) + fvc::div(phi, K)
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
        - fvm::laplacian(thermo.alpha()/rho, hs)
        //+ mstm().multiSpeciesHeatSource()
        ==
          chemistrySh
        + fvOptions(rho, hs)
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
