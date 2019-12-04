% Script to evaluate the particle contact with the stl-walls:
% 1. direct contacts of particles with the walls
% 2. indirect contacts of particles through contact chains

% Clear workspace and command screen
clear;
clc;
close all;
% Contact information is read from the custom and local dump files
WallContactData = readdump_wallgranlocal('../DEM/post/dump.wallgranlocal_restart');
PairContactData = readdump_pairgranlocal('../DEM/post/dump.pairgranlocal_restart');
AtomData = readdump_all('../DEM/post/dump.liggghts_restart');

Ntimesteps = length(AtomData.timestep);      % Number of timesteps [-]
WallContactFlag = zeros(AtomData.Natoms(1),2,Ntimesteps); % Flag array indicating wall contact for each particle id
ParticleRadius = zeros(AtomData.Natoms(1),2,Ntimesteps); % Particle radius array

TimeStepSize = 1e-6;     % [s]
Time = AtomData.timestep*TimeStepSize;

% Particle surface area at each timestep of all particles 
% with direct or indirect wall contact [m²]
ActiveParticleSurface = zeros(Ntimesteps,1);
% Count contacts
ParticleContacts = zeros(Ntimesteps,1);

for i = 1 : Ntimesteps

	% Set particle radius array 
	for j = 1 : AtomData.Natoms(i)
		id = AtomData.atom_data(j,1,i);
		ParticleRadius(id,1,i) = id;
		ParticleRadius(id,2,i) = AtomData.atom_data(j,19,i)*0.01;
	endfor

	% Find direct wall contacts
	for j = 1 : WallContactData.Nwallcontacts(i)
		% Get atom id of direct wall contact
		id = WallContactData.wallcontact_data(j,10,i);
		WallContactFlag(id,1,i) = id;
		WallContactFlag(id,2,i) = 1;
		% Add direct wall contact to active surface area
		ActiveParticleSurface(i) = ActiveParticleSurface(i) ...
                                  + ParticleRadius(id,2,i)^2*4.0*pi;
		ParticleContacts(i) = ParticleContacts(i) + 1;
	endfor
	
    LoopContactCount = 1;
	% Find indirect contacts
	while(LoopContactCount > 0)
		LoopContactCount = 0;
		% Check pair contact list for additonal indirect contacts
		for j = 1 : PairContactData.Npaircontacts(i)
			id1 = PairContactData.paircontact_data(j,2,i);
			id2 = PairContactData.paircontact_data(j,3,i);
			if (WallContactFlag(id1,2,i) == 0 && WallContactFlag(id2,2,i) == 1)
				WallContactFlag(id1,1,i) = id1;
				WallContactFlag(id1,2,i) = 1;
				LoopContactCount = LoopContactCount + 1;
				% Add indirect wall contact to active surface area
				ActiveParticleSurface(i) = ActiveParticleSurface(i) ...
										  + ParticleRadius(id1,2,i)^2*4.0*pi;
				ParticleContacts(i) = ParticleContacts(i) + 1;
			elseif (WallContactFlag(id1,2,i) == 1 && WallContactFlag(id2,2,i) == 0)
				WallContactFlag(id2,1,i) = id2;
				WallContactFlag(id2,2,i) = 1;
				LoopContactCount = LoopContactCount + 1;
				% Add indirect wall contact to active surface area
				ActiveParticleSurface(i) = ActiveParticleSurface(i) ...
										  + ParticleRadius(id2,2,i)^2*4.0*pi;
				ParticleContacts(i) = ParticleContacts(i) + 1;
			endif
		endfor
	endwhile
endfor
			
fid = fopen ("particleContacts.txt", "w");
fprintf(fid, "File contains the active particle surface and all particle contacts with the collector surface over the simulated time\n");
fprintf(fid, "Time [s]	Active Particle Surface [m²]	# of Contacting Particles [-]	Time Step [-]\n");
for i = 1 : Ntimesteps
	fprintf(fid, "%12.3e %12.3e %10i %10i\n", Time(i), ActiveParticleSurface(i), ParticleContacts(i), AtomData.timestep(i));
endfor




