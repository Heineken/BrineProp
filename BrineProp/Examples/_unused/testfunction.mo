within BrineProp.Examples._unused;
function testfunction
 input Real T;
 output Real h;
algorithm
  h := ln(T);
  Modelica.Utilities.Streams.print(String(T));
end testfunction;
