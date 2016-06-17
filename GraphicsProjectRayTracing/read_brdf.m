function [B,format]= read_brdf(filename, channel)
% [B, format] = read_brdf(filename, channel)
%
% Output: B - BRDF ( num_channel x dim0 x dim1 x dim2 x dim3 ) , format - struct with tabular format parameters
% Input: filename, channel ( 0 - 'red' , 1 - 'green', 2 - 'blue', -1 - 'All')
%
% Written by Addy Ngan - MIT
%
% Ver 1.01 7/20/2005
% Ver 1.0 6/19/2005

fid = fopen(filename, 'rb');
format.dims = fread(fid,4,'uint32');
% Read in dimensions of the tabulated BRDF

% Standard Parameterization: Theta_In x Theta_Out x Phi_Diff x Phi_In
% Theta_In and Theta_Out range from 0 to pi/2, Phi_Diff range from 0 to pi
% (if half_data == 0) or 2pi (if half_data == 1), Phi_In from 0 to 2pi
% for anisotropic materials.

[modes, count1] = fread(fid,10,'int32');
[expon, count2] = fread(fid,1,'double');

if (count1~=10 || count2~=1)
    error('File error');
else
    format.paramType = modes(3);  % 0 - Rusinkiewicz Parameterization, 1 - Standard Parameterization
    format.binType = modes(4); % 0 - linear binning interval, 1 - reserved
    format.half_data = modes(6)>0; % 0 - Phi_Diff range [0 2pi] , 1 - Phi_Diff range [0 pi]
    format.num_channels = modes(7); % Always 3 for our data (R,G,B)
    format.expon = expon; % reserved
end

num_channels = format.num_channels;

if channel~=-1
    fseek(fid, channel*prod(format.dims)*4, 0);
    num_channels=1;
end

B = fread(fid,num_channels*prod(format.dims),'float=>single');
fclose(fid);

if prod(size(B))<num_channels*prod(format.dims)
    B=0;
    error('Loading failed');
end

% Reshaping the data into desired dimensions
B=reshape(B,[fliplr(format.dims') num_channels]);
p = length(size(B)):-1:1;
B = permute(B,p);
if num_channels ==1
    B = reshape(B, [1 size(B)]);
end
    



