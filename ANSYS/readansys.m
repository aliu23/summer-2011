function [ var ] = readansys( file )

f=textread(file);

switch f(1)
  
  case 1
    
      N=f(2); %number of steps
      B=8;    %number of lines until next component
      angles=linspace((f(3).*180/pi),f(4).*180/pi,N);   %creates a vector of angles
      Fxmx=f(5:B:end);     %Force x (maxwell stress tensor)
      Fxvw=f(6:B:end);     %Force x (virtual work)
      Fymx=f(7:B:end);     %Force y (maxwell stress tensor)
      Fyvw=f(8:B:end);     %Force y (virtual work)
      Tzmx=f(11:B:end);    %Torque z local (maxwell stress tensor)
      Tzvw=f(12:B:end);    %Torque z local (virtual work)
      
      var.description = '2D torques and forces vs angles';
      var.angles = angles;
      var.force_mx = [Fxmx Fymx];
      var.force_vw = [Fxvw Fyvw];
      var.torque_mx = Tzmx;
      var.torque_vw = Tzvw;
   
  case 2
      
      N=f(2); %number of steps
      B=6;    %number of lines until next component
      distancex=linspace(f(3),f(4),N);   %creates a vector of distances
      Fxmx=f(5:B:end);     %Force x (maxwell stress tensor)
      Fxvw=f(6:B:end);     %Force x (virtual work)
      Fymx=f(7:B:end);     %Force y (maxwell stress tensor)
      Fyvw=f(8:B:end);     %Force y (virtual work)
      Tzmx=f(9:B:end);     %Torque z local (maxwell stress tensor)
      Tzvw=f(10:B:end);    %Torque z local (virtual work)
      
      var.description = '2D torques and forces vs displacement';
      var.distance = distancex;
      var.force_mx = [Fxmx Fymx];
      var.force_vw = [Fxvw Fyvw];
      var.torque_mx = Tzmx;
      var.torque_vw = Tzvw;

  case 3
      
      N=f(2); %number of steps
      B=8;    %number of lines until next component
      distancex=linspace(f(3),f(4),N);   %creates a vector of distances
      Fxmx=f(5:B:end);     
      Fxvw=f(6:B:end);    
      Fymx=f(7:B:end);     
      Fyvw=f(8:B:end);
      Fzmx=f(9:B:end);
      Fzvw=f(10:B:end);
      Tzmx=f(11:B:end);   
      Tzvw=f(12:B:end);  
      
      var.description = '3D torques and forces vs displacement';
      var.distance = distancex;
      var.force_mx = [Fxmx Fymx Fzmx];
      var.force_vw = [Fxvw Fyvw Fzvw];
      var.torque_mx = Tzmx;
      var.torque_vw = Tzvw;

  otherwise
    
    error('Unknown data type')

end



end

