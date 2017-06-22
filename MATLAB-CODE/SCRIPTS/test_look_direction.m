for direction = [1,2,3];
    fname=sprintf('test_look_dir_%d',direction);
    back_project_synthetics(15,1,5,5,logspace(-4,-1,3),[2],500,direction,fname);
end