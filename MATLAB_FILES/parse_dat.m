% Copyright (c) 2019 TainiTec Ltd.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [d,t] = parse_dat(dat_fn)
  
%  printf('\nLoading DAT file: %s\n\n', dat_fn);
  
  % Parameters (TODO load from configuration file)
  N_CHANNELS = 16; % number of channels per transmitter
  SAMPLE_RATE = 250.40; % Hz
  %SAMPLE_RATE = 19531.25; % Hz
  MVOLT_RANGE = 13.0; % mV
  INTEGER_RANGE = 2^12-1; % the maximum integer value from TAINI
  MISSED_PACKET_SYMBOL = 2^15-1; % the value substituted into the DAT file during periods of packet loss
  
  % Load the raw data as a single dimensional array
  dat_fid = fopen(dat_fn);
  raw_dat = fread(dat_fid,'int16');
  
  % Calculate the signal length
  n_samples = length(raw_dat) / N_CHANNELS;
  
  % De-multiplex the channels
  sample_integer_values = reshape(raw_dat, N_CHANNELS, n_samples)';
  
  % Convert the integer values to milli-volts
  d = sample_integer_values / INTEGER_RANGE * MVOLT_RANGE;
  
  % Form the time vector
  t = (0:(n_samples-1)) / SAMPLE_RATE;
  
  % Remove any missed packet samples (replace with nan)
  d(sample_integer_values == MISSED_PACKET_SYMBOL) = 0;
  
end
