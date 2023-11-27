%LOAD_VARY_QUANT_LEVELS_DATA
function [scfile, measidx, p, rhos, ges_suffix, imindex] = load_vary_quant_levels_data(quant_level)

switch quant_level
    case 4
        scfile     = 2028;
        measidx    = [1];
        p          = 2;
        rhos       = [1.5]*1e2;
        ges_suffix = '_touch_4bit';
        imindex    = 'e';
    case 4.4
        scfile     = 2020;
        measidx    = [2];
        p          = 3;
        rhos       = [6]*1e2;
        ges_suffix = '_touch_4d4bit';
        imindex    = 'd';
    case 4.9
        scfile     = 2006;
        measidx    = [1];
        p          = 3;
        rhos       = [4]*1e2;
        ges_suffix = '_touch_4d9bit';
        imindex    = 'c';
    case 5.2
        scfile     = 2032;
        measidx    = [4];
        p          = 3;
        rhos       = [6]*1e2;
        ges_suffix = '_touch_5d2bit';
        imindex    = 'b';
    case 5.6
        scfile     = 2007;
        measidx    = [10];
        p          = 3;
        rhos       = [6]*1e2;
        ges_suffix = '_touch_5d6bit';
        imindex    = 'a';
    otherwise 
        disp('unsupported quantization level.');
end
