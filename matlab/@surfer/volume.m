function v = volume(obj)
% VOLUME of surfer object ( 1/3 * \int_S r \cdot n \, dS ) 
%

v = (sum(obj.r.*obj.n,1)*obj.wts)/3;
