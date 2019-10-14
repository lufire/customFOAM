function [varargout] = readdump_pairgranlocal(varargin)
% Reads all timesteps from a LAMMPS dump file. 
% Input is dump file name with path
% Output is in the form of a structure with following variables
% .timestep     --> Vector containing all time steps
% .Natoms       --> Vector containing number of atoms at each time step
% .x_bound      --> [t,2] array with xlo,xhi at each time step
% .y_bound      --> [t,2] array with ylo,yhi at each time step
% .z_bound      --> [t,2] array with zlo,zhi at each time step
% .atom_data    --> 3 dimensional array with data at each time step stored
%                   as atomdata(:,:,t)
% Example
%       data = readdump_all('dump.LAMMPS'); 
%
% See also readdump_one, scandump
%
%  Author :  Arun K. Subramaniyan
%            sarunkarthi@gmail.com
%            http://web.ics.purdue.edu/~asubrama/pages/Research_Main.htm
%            School of Aeronautics and Astronautics
%            Purdue University, West Lafayette, IN - 47907, USA.

try
    dump = fopen(varargin{1},'r');
catch
    error('Dumpfile not found!');
end

i=1;
while feof(dump) == 0
    id = fgetl(dump);
	if (strcmpi(id,'ITEM: TIMESTEP'))
        timestep(i) = str2num(fgetl(dump));
    else
		if (strcmpi(id,'ITEM: NUMBER OF ENTRIES'))
			Npaircontacts(i) = str2num(fgetl(dump));
		else
			if (strcmpi(id(1:16),'ITEM: BOX BOUNDS'))
				x_bound(i,:) = str2num(fgetl(dump));
				y_bound(i,:) = str2num(fgetl(dump));
				z_bound(i,:) = str2num(fgetl(dump));
			else
				if (strcmpi(id(1:13),'ITEM: ENTRIES'))
					for j = 1 : 1: Npaircontacts(i)
						paircontact_data(j,:,i) = str2num(fgetl(dump));
					end
					i=i+1;
				end
			end 
		end
	end
end
%----------Outputs-------------
%OUTPUTS IN SAME VARIABLE STRUCTURE
varargout{1}.timestep = timestep;
varargout{1}.Npaircontacts = Npaircontacts;
varargout{1}.x_bound = x_bound;
varargout{1}.y_bound = y_bound;
varargout{1}.z_bound = z_bound;
varargout{1}.paircontact_data = paircontact_data;
%------------------------------
fclose(dump);
