clear all;
close all;

dir = './0';

% ----- start Calculating ----------

dir_X=sprintf([dir '/ccx']);
dir_Y=sprintf([dir '/ccy']);
dir_Z=sprintf([dir '/ccz']);
dirv=sprintf([dir '/V']);

delimiterIn=' ';
headerlinesIn=22;

X = importdata(dir_X,delimiterIn,headerlinesIn);
Y = importdata(dir_Y,delimiterIn,headerlinesIn);
Z = importdata(dir_Z,delimiterIn,headerlinesIn);
v = importdata(dirv,delimiterIn,headerlinesIn);

dir_alpha=sprintf([dir '/alpha.liq']);
dir_alpha_field=sprintf([dir '/alpha.liq']);

alpha = importdata(dir_alpha,delimiterIn,headerlinesIn);

[n_cell,~] = size(X.data);

delta = 0.3;

overwrite = true;

r = 20e-6;
x = 0;
y = 0;
z = 0;

delta = delta*max(power(v.data(:),1/3));

for i=1:1:n_cell       

    r_xyz = sqrt(power(X.data(i)-x,2) + power(Y.data(i)-y,2) + power(Z.data(i)-z,2));

    var = (r_xyz - r)/(delta + 1e-16);
  
    if var > 1

        if overwrite == true

            alpha.data(i) = 0;

        end;

    elseif var < -1 

            alpha.data(i) = 1;

    else
         if overwrite == false && alpha.data(i) ~= 1

            alpha.data(i) = -0.5*var + 0.5;

         end;
         
         if overwrite == true

            alpha.data(i) = -0.5*var + 0.5;

         end;
         
    end;
    
end;

fileID = fopen(dir_alpha_field,'w');

for i=1:1:9    
    fprintf(fileID,'%s\n',X.textdata{i});
end
fprintf(fileID,'%s\n','    version     2.3;');
fprintf(fileID,'%s\n','    format      ascii;');
fprintf(fileID,'%s\n','    class       volScalarField;');
fprintf(fileID,'%s\n','    location    "0";');
fprintf(fileID,'%s\n','    object      alpha.liq;');
fprintf(fileID,'%s\n','}');
fprintf(fileID,'%s\n','// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //');

fprintf(fileID,'\n%s\n\n%s\n%i\n%s\n','dimensions      [0 0 0 0 0 0 0];','internalField   nonuniform List<scalar>',n_cell,'(');

for i=1:1:n_cell
    fprintf(fileID,'%f\n',alpha.data(i));
end

fprintf(fileID,'%c\n',');');
fprintf(fileID,'%s\n','boundaryField');
fprintf(fileID,'%c\n','{');
fprintf(fileID,'%s\n','    other');
fprintf(fileID,'%s\n','    {');
fprintf(fileID,'%s\n','        type            zeroGradient;');
fprintf(fileID,'%s\n','    }');
fprintf(fileID,'%s\n','    symmetryPlane');
fprintf(fileID,'%s\n','    {');
fprintf(fileID,'%s\n','        type            symmetry;');
fprintf(fileID,'%s\n','    }');
fprintf(fileID,'%s\n\n','}');
fprintf(fileID,'%s','// ************************************************************************* //');

fclose(fileID);

