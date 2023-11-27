%LOAD_GESTURE_DATA
function [scfile, measidx, rhos, ges_suffix] = load_gesture_data(gesture)

switch lower(gesture)
    case 'slide'
        scfile     = 1912;
        measidx    = [21,27,28,32,41,42];
        rhos       = [6,4,5,6,6,6]*1e2;
        ges_suffix = '_slide1_6p';
    case 'scroll'
        scfile     = 1925;
        measidx    = [9,11,13,16,22,23];
        rhos       = [5,6,6,13,9,11]*1e2;
        ges_suffix = '_scroll2_6p';
    case 'pinch'
        scfile     = 1924;
        measidx    = [21,26,27,28,29,30];
        rhos       = [5,5,4,5,5,5]*1e2;
        ges_suffix = '_pinch3_6p';
    case 'swipe'
        scfile     = 1916;
        measidx    = [19,20,21,24,26,31];
        rhos       = [6,5,5,10,8,11]*1e2;
        ges_suffix = '_swipe4_6p';
    case 'rotate'
        scfile     = 1927;
        measidx    = [54,59,62,67,68,74];
        rhos       = [6,5,10,6,6,5]*1e2;
        ges_suffix = '_rotate5_6p';
    otherwise 
        disp('unsupported gesture.');
end
