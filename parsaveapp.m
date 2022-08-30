%% Function for save mat file in parallel for loop
function parsaveapp(fname,x)
    save(fname,'x','-append');
end