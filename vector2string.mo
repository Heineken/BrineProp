within Brine;
function vector2string "create string from vector"
  input Real[:] vec;
  input Boolean neat=true "Align and crop elements";
  output String str;
algorithm
  str:=if neat then "\t{" else "{";
  for i in 1: size(vec,1) loop
    str:=str + (if i>1 then "," else "");
    if neat then
      str:=str + String(vec[i],false,12,significantDigits=5);
    else
      str:=str + String(vec[i]);
    end if;
  end for;
  str:=str + "}";
end vector2string;
