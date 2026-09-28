function key = encode_state(x, base)
%ENCODE_STATE Encode nonnegative integer vector as a uint64 scalar key.
x = uint64(x(:).'); base = uint64(base);
key = uint64(0); mult = uint64(1);
for k = 1:numel(x)
    key = key + x(k)*mult;
    mult = mult*base;
end
end
