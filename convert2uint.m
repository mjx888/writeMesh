function [convertedArray, selectedType] = convert2uint( inputArray )
% convert2uint: Converts double array to the smallest possible 
% unsigned integer type.
%
% [convertedArray, selectedType] = smart_uint_converter(inputArray)
% 
% Inputs:
%   inputArray - An array of doubles (can be any dimension).
% 
% Outputs:
%   convertedArray - The array converted to uint8, uint16, uint32, or 
%                    uint64.
%   selectedType   - A string indicating which type was selected.
% 
% Logic:
%   - If max value <= 255, converts to uint8.
%   - If max value <= 65535, converts to uint16.
%   - If max value <= 4294967295, converts to uint32.
%   - Otherwise, converts to uint64.
% 
% Notes:
%   - Negative values in the input will be converted to 0.
%   - Non-integer values will be rounded to the nearest integer.
%   - Values larger than intmax('uint64') will saturate to that maximum.
%

    % Validate input
    if ~isnumeric(inputArray)
        error('Input must be a numeric array.');
    end

    % Check if input is already an unsigned integer type
    if isa(inputArray, 'uint8') || isa(inputArray, 'uint16') || ...
       isa(inputArray, 'uint32') || isa(inputArray, 'uint64')
        convertedArray = inputArray;
        selectedType = class(inputArray);
        return;
    end

    % Handle empty array case
    if isempty(inputArray)
        convertedArray = uint8([]);
        selectedType = 'uint8';
        return;
    end

    % Find the maximum value in the array ignoring NaNs
    % We use double() on intmax to ensure safe comparison
    maxVal = max(inputArray(:));

    % Determine appropriate type based on the maximum value
    if maxVal <= double(intmax('uint8'))
        convertedArray = uint8(inputArray);
        selectedType = 'uint8';
        
    elseif maxVal <= double(intmax('uint16'))
        convertedArray = uint16(inputArray);
        selectedType = 'uint16';
        
    elseif maxVal <= double(intmax('uint32'))
        convertedArray = uint32(inputArray);
        selectedType = 'uint32';
        
    else
        % If it exceeds uint32, default to the largest available: uint64
        convertedArray = uint64(inputArray);
        selectedType = 'uint64';
        
        % Optional check: If the value is actually larger than uint64 can hold
        if maxVal > 1.844674407370955e+19
            warning('Values in array exceed the maximum limit of uint64. Data saturated.');
        end
    end
end