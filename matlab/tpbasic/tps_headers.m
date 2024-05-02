function s = tps_headers(f)
% function s = tps_headers(f)
%---
% reads headers from CFD files or XML files
% and transcripts into a more readable structure

s = tps_read(f,'header');