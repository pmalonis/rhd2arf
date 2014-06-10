function arfwriteatt(filename,location,attname,attvalue)
    if numel(attvalue) > 1 || ~isnumeric(attvalue)
        h5writeatt(filename,location,attname,attvalue)        
    else
        if isinteger(attvalue)
            type_id = H5T.copy('H5T_NATIVE_INT');
        else
            type_id = H5T.copy('H5T_NATIVE_DOUBLE');
        end
        fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
        acpl_id = H5P.create('H5P_ATTRIBUTE_CREATE');
        loc_id = H5O.open(fid,location,'H5P_DEFAULT');
        space_id = H5S.create('H5S_SCALAR');
        attr_id = H5A.create(loc_id,attname,type_id,space_id,acpl_id);
        H5A.write(attr_id,'H5ML_DEFAULT',attvalue)
    end
end