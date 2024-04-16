#include "H5Writer.h"

// Constructor
H5Writer::H5Writer(std::string fileName)
{
    // Create a new file using default properties
    file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Initialize the dataset dimensions
    dims_1D[0] = 0;
    dims_2D[0] = 0;
    dims_2D[1] = 0;
}

// Destructor
H5Writer::~H5Writer()
{
    // Close the file
    H5Fclose(file_id);
}

// Create a new group with method chaining
H5Writer& H5Writer::createGroup(std::string groupName)
{
    // Create a new group
    group_id = H5Gcreate(file_id, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Return a reference to the current object
    return *this;
}

// Close the current group, with method chaining
H5Writer& H5Writer::closeGroup()
{
    // Close the current group
    H5Gclose(group_id);

    // Return a reference to the current object
    return *this;
}

// Write a H5T_STD_I64LE dataset, with method chaining
H5Writer& H5Writer::writeDataset(std::string datasetName, uint64_t* data, int &rank, hsize_t* dims)
{
    // Create a new dataspace
    dataspace_id = H5Screate_simple(rank, dims, NULL);

    // Create a new dataset
    dataset_id = H5Dcreate(file_id, datasetName.c_str(), H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write the data to the dataset
    status = H5Dwrite(dataset_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    // Close the dataset
    H5Dclose(dataset_id);

    // Close the dataspace
    H5Sclose(dataspace_id);

    // Return a reference to the current object
    return *this;
}

// Write a H5T_IEEE_F64LE dataset, with method chaining
H5Writer& H5Writer::writeDataset(std::string datasetName, double* data, int &rank, hsize_t* dims)
{
    // Create a new dataspace
    dataspace_id = H5Screate_simple(rank, dims, NULL);

    // Create a new dataset
    dataset_id = H5Dcreate(file_id, datasetName.c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write the data to the dataset
    status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    // Close the dataset
    H5Dclose(dataset_id);

    // Close the dataspace
    H5Sclose(dataspace_id);

    // Return a reference to the current object
    return *this;
}