#ifndef H5WRITER_H_
#define H5WRITER_H_

#include <hdf5.h>
#include <hdf5_hl.h>
#include <cstdint>
#include <string>
#include <pcgnslib.h>

class H5Writer
{
    private:
        hid_t file_id;
        herr_t status;
        hsize_t dims_1D[1];
        hsize_t dims_2D[2];
        hid_t group_id;
        hid_t dataspace_id;
        hid_t dataset_id;
        std::string groupName;
        std::string datasetName;

        // Create  a new group with method chaining
        H5Writer& createGroup(std::string groupName);
        // Close the current group, with method chaining
        H5Writer& closeGroup();
        // Write a H5T_STD_I64LE dataset, with method chaining
        H5Writer& writeDataset_INT(std::string datasetName, int64_t* data, int &rank, hsize_t* dims);
        // Write a H5T_IEEE_F64LE dataset, with method chaining
        H5Writer& writeDataset_REAL(std::string datasetName, double* data, int &rank, hsize_t* dims);
    public:
        // Constructor
        H5Writer(std::string fileName);
        // Destructor
        ~H5Writer();
        // Write the SOD2D file
        void writeSOD2D( uint64_t &nElem, uint64_t &nBelem, uint64_t &nPoin, int &nNode, int &pOrder,
                         cgsize_t** &connecSOD2D, uint64_t** &connecBoundSOD2D, double* &coordSOD2D );
};

#endif //! H5WRITER_H_