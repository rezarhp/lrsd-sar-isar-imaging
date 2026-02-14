function SAR_cmp = sar_cmap
% SAR_cmp = sar_cmap
%
% INPUTS: None
%
% OUTPUTS
%
% SAR_cmp : A colormap for viewing SAR imagery
%
% This routine generates a colormap, useful for displaying SAR imagery. The
% colormap can subsequently be invoked by the MATLAB command
%
% colormap(SAR_cmp)

% W. C. Karl 6/11/97 based on code from W. W. Irving.

SAR_cmp = zeros(256,3);

% Red
SAR_cmp(80:191,1) = 0.9961 * [0:111]'/111;
SAR_cmp(192:256,1) = 0.9961 * ones(65,1);

% Green
SAR_cmp(3:256,2) = 0.9844 * [0:253]'/253;
 
% Blue
SAR_cmp(1:20,3) = 0.2344 * [0:19]'/19;
SAR_cmp(21:30,3) = 0.2344 + 0.039*[0:9]'/9;
SAR_cmp(31:40,3) = 0.2734 - 0.039*[0:9]'/9;
SAR_cmp(41:80,3) = 0.2344 - 0.2344*[0:39]'/39;
SAR_cmp(180:245,3) = 0.9961 * [0:65]'/65;
SAR_cmp(246:256,3) = 0.9961 * ones(11,1);