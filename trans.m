function [mat] = trans(x,y,z,a,b,g)
    mat = eye(4);
    mat(1:3,1:3) = eul2rotm([a,b,g],'XYZ');
    mat(1:3,4) = [x;y;z];
end

