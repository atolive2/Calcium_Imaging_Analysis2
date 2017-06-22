%change_size code --------------------------------------------------
%copy paste this code into a file called change_size.m
function  change_size(objHandel, evt, annotation_handle, orig_pos)
    slider_value = get(objHandel,'Value');
    new_pos = orig_pos;
    new_pos(3:4) = orig_pos(3:4)*slider_value;
    set(annotation_handle, 'position', new_pos)