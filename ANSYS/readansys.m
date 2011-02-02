function [ var ] = readansys( file )

f=textread(file);

switch f(1)
  
  case 1
    
   N=f(2); %number of steps
   B=8;    %number of lines until next component
   angles=linspace((f(3).*180/pi),f(4).*180/pi,N);   %creates a vector of angles 
   Fxmx=f(5:B:end);     %Force x (maxwell stress tensor)
   Fxvw=f(6:B:end);     %Force x (virtual work)
   Fymx=f(7:B:end);
   Fyvw=f(8:B:end);
   Tzmx=f(11:B:end);    %Torque z local (maxwell stress tensor)
   Tzvw=f(12:B:end);    %Torque z local (virtual work)
   
   var.description = '2D torques vs angles';
   var.angles = angles;
   var.force_mx = [Fxmx Fymx];
   var.force_vw = [Fxvw Fyvw];
   var.torque_mx = Tzmx;
   var.torque_vw = Tzvw;
   
  case 2
    
  case 3
    
  otherwise
    
    error('Unknown data type')

end



end

