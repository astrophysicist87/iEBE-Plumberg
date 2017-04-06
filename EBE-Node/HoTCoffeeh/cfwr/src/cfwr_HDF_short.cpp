#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>
#include<algorithm>

#include "H5Cpp.h"

#include "cfwr.h"
#include "lib.h"

using namespace std;

const int RANK2D = 2;

//*******************************************
// HDF array for resonances
//*******************************************
/////////////////////////////////////////////
int CorrelationFunction::Initialize_resonance_HDF_array()
{
const int chunk_size = n_pT_pts * n_pphi_pts * /*n_pY_pts **/ qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

	double * resonance_chunk = new double [chunk_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << path << "/resonance_spectra.h5";
	H5std_string RESONANCE_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string RESONANCE_DATASET_NAME("ra");

	//bool file_does_not_already_exist = !fexists(filename_stream_ra.str().c_str());
	bool file_does_not_already_exist = true;	//force full initialization for timebeing...

	try
    {
		Exception::dontPrint();
	
		resonance_file = new H5::H5File(RESONANCE_FILE_NAME, H5F_ACC_TRUNC);

		DSetCreatPropList cparms;
		hsize_t chunk_dims[RANK2D] = {1, chunk_size};
		cparms.setChunk( RANK2D, chunk_dims );

		hsize_t dims[RANK2D] = {NchosenParticle + 1, chunk_size};
		resonance_dataspace = new H5::DataSpace (RANK2D, dims);

		resonance_dataset = new H5::DataSet( resonance_file->createDataSet(RESONANCE_DATASET_NAME, PredType::NATIVE_DOUBLE, *resonance_dataspace, cparms) );

		hsize_t count[RANK2D] = {1, chunk_size};
		hsize_t dimsm[RANK2D] = {1, chunk_size};
		hsize_t offset[RANK2D] = {0, 0};

		resonance_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
		if (file_does_not_already_exist)
		{
			*global_out_stream_ptr << "HDF resonance file doesn't exist!  Initializing to zero..." << endl;

			for (int ir = 0; ir < NchosenParticle + 1; ++ir)
			for (int ipY = 0; ipY < n_pY_pts; ++ipY)
			{
				//offset[0] = ir;
				offset[0] = ir * n_pY_pts + ipY;
				resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
	
				for (int iidx = 0; iidx < chunk_size; ++iidx)
					resonance_chunk[iidx] = 0.0;
	
				//initialize everything with zeros
				resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
			}
		}
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] resonance_chunk;

	return (0);
}

int CorrelationFunction::Initialize_target_thermal_HDF_array()
{
const int chunk_size = n_pT_pts * n_pphi_pts * /*n_pY_pts **/ qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

	double * tta_chunk = new double [chunk_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << path << "/target_thermal_moments.h5";
	H5std_string TARGET_THERMAL_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string TARGET_THERMAL_DATASET_NAME("tta");

	//bool file_does_not_already_exist = !fexists(filename_stream_ra.str().c_str());
	bool file_does_not_already_exist = true;	//force full initialization for timebeing...

	try
    {
		Exception::dontPrint();
	
		tta_file = new H5::H5File(TARGET_THERMAL_FILE_NAME, H5F_ACC_TRUNC);

		DSetCreatPropList cparms;
		hsize_t chunk_dims[1] = {chunk_size};
		cparms.setChunk( 1, chunk_dims );

		hsize_t dims[1] = {chunk_size};
		tta_dataspace = new H5::DataSpace (1, dims);

		tta_dataset = new H5::DataSet( tta_file->createDataSet(TARGET_THERMAL_DATASET_NAME, PredType::NATIVE_DOUBLE, *tta_dataspace, cparms) );

		hsize_t count[1] = {chunk_size};
		hsize_t dimsm[1] = {chunk_size};
		hsize_t offset[1] = {0};

		tta_memspace = new H5::DataSpace (1, dimsm, NULL);
		if (file_does_not_already_exist)
		{
			*global_out_stream_ptr << "HDF thermal target moments file doesn't exist!  Initializing to zero..." << endl;

			for (int ipY = 0; ipY < n_pY_pts; ++ipY)
			{

				offset[0] = ipY;
				tta_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
	
				for (int iidx = 0; iidx < chunk_size; ++iidx)
					tta_chunk[iidx] = 0.0;
	
				//initialize everything with zeros
				tta_dataset->write(tta_chunk, PredType::NATIVE_DOUBLE, *tta_memspace, *tta_dataspace);
			}
		}
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] tta_chunk;

	return (0);
}

/////////////////////////////////////////////

int CorrelationFunction::Open_resonance_HDF_array(string resonance_local_file_name)
{
const int chunk_size = n_pT_pts * n_pphi_pts * /*n_pY_pts **/ qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

	ostringstream filename_stream_ra;
	filename_stream_ra << path << "/" << resonance_local_file_name;
	H5std_string RESONANCE_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string RESONANCE_DATASET_NAME("ra");

	try
    {
		Exception::dontPrint();
	
		hsize_t dimsm[RANK2D] = {1, chunk_size};
		hsize_t dims[RANK2D] = { (NchosenParticle + 1) * n_pY_pts, chunk_size};

		resonance_dataspace = new H5::DataSpace (RANK2D, dims);
		resonance_file = new H5::H5File(RESONANCE_FILE_NAME, H5F_ACC_RDWR);
		resonance_dataset = new H5::DataSet( resonance_file->openDataSet( RESONANCE_DATASET_NAME ) );
		resonance_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }


	return (0);
}

/////////////////////////////////////////////

int CorrelationFunction::Open_target_thermal_HDF_array()
{
const int chunk_size = n_pT_pts * n_pphi_pts * /*n_pY_pts **/ qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

	ostringstream filename_stream_ra;
	filename_stream_ra << path << "/target_thermal_moments.h5";
	H5std_string TARGET_THERMAL_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string TARGET_THERMAL_DATASET_NAME("tta");

	try
    {
		Exception::dontPrint();
	
		hsize_t dimsm[1] = {chunk_size};
		hsize_t dims[1] = {chunk_size * n_pY_pts};

		tta_dataspace = new H5::DataSpace (1, dims);
		tta_file = new H5::H5File(TARGET_THERMAL_FILE_NAME, H5F_ACC_RDWR);
		tta_dataset = new H5::DataSet( tta_file->openDataSet( TARGET_THERMAL_DATASET_NAME ) );
		tta_memspace = new H5::DataSpace (1, dimsm, NULL);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	return (0);
}

/////////////////////////////////////////////

int CorrelationFunction::Close_resonance_HDF_array()
{
	try
    {
		Exception::dontPrint();
	
		resonance_memspace->close();
		resonance_dataset->close();
		resonance_file->close();
		delete resonance_memspace;
		delete resonance_file;
		delete resonance_dataset;
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	return (0);
}

/////////////////////////////////////////////

int CorrelationFunction::Close_target_thermal_HDF_array()
{
	try
    {
		Exception::dontPrint();
	
		tta_memspace->close();
		tta_dataset->close();
		tta_file->close();
		delete tta_memspace;
		delete tta_file;
		delete tta_dataset;
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	return (0);
}

/////////////////////////////////////////////

int CorrelationFunction::Set_resonance_in_HDF_array(int local_pid, int ipY, double * resonance_array_to_use)
{
const int chunk_size = n_pT_pts * n_pphi_pts * /*n_pY_pts **/ qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

	double * resonance_chunk = new double [chunk_size];

	int local_icr;
	if (local_pid == target_particle_id)
		local_icr = 0;
	else
		local_icr = lookup_resonance_idx_from_particle_id(local_pid) + 1;	//note shift

	try
    {
		Exception::dontPrint();
		hsize_t offset[RANK2D] = {local_icr * n_pY_pts + ipY, 0};
		hsize_t count[RANK2D] = {1, chunk_size};				// == chunk_dims
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		resonance_dataset->write(resonance_array_to_use, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] resonance_chunk;

	return (0);
}

/////////////////////////////////////////////

int CorrelationFunction::Get_resonance_from_HDF_array(int local_pid, int ipY, double * resonance_array_to_fill)
{
const int chunk_size = n_pT_pts * n_pphi_pts * /*n_pY_pts **/ qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

	double * resonance_chunk = new double [chunk_size];

	//store target particle moments in first slot of resonance array
	int local_icr;
	if (local_pid == target_particle_id)
		local_icr = 0;
	else
		local_icr = lookup_resonance_idx_from_particle_id(local_pid) + 1;	//note shift

	try
    {
		Exception::dontPrint();
		hsize_t offset[RANK2D] = {local_icr * n_pY_pts + ipY, 0};
		hsize_t count[RANK2D] = {1, chunk_size};				// == chunk_dims
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		resonance_dataset->read(resonance_array_to_fill, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] resonance_chunk;

	return (0);
}

/////////////////////////////////////////////

int CorrelationFunction::Set_target_thermal_in_HDF_array(int ipY, double * tta_array_to_use)
{
const int chunk_size = n_pT_pts * n_pphi_pts * /*n_pY_pts **/ qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

	double * tta_chunk = new double [chunk_size];

	try
    {
		Exception::dontPrint();
	
		hsize_t offset[1] = {ipY};
		hsize_t count[1] = {chunk_size};				// == chunk_dims
		tta_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		tta_dataset->write(tta_array_to_use, PredType::NATIVE_DOUBLE, *tta_memspace, *tta_dataspace);
   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] tta_chunk;

	return (0);
}


/////////////////////////////////////////////

//int CorrelationFunction::Get_target_thermal_from_HDF_array(double ******* tta_to_fill)
int CorrelationFunction::Get_target_thermal_from_HDF_array(int ipY, double * tta_to_fill)
{
const int chunk_size = n_pT_pts * n_pphi_pts * /*n_pY_pts **/ qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

	double * tta_chunk = new double [chunk_size];

	try
    {
		Exception::dontPrint();
		hsize_t offset[1] = {ipY};
		hsize_t count[1] = {chunk_size};				// == chunk_dims
		tta_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		tta_dataset->read(tta_to_fill, PredType::NATIVE_DOUBLE, *tta_memspace, *tta_dataspace);
   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] tta_chunk;

	return (0);
}


/////////////////////////////////////////////

int CorrelationFunction::Copy_chunk(int current_resonance_index, int reso_idx_to_be_copied)
{
const int chunk_size = n_pT_pts * n_pphi_pts * /*n_pY_pts **/ qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

	double * resonance_chunk = new double [chunk_size];

	int local_icr1, local_icr2;
	if (current_resonance_index == target_particle_id)
		local_icr1 = 0;
	else
		local_icr1 = lookup_resonance_idx_from_particle_id(current_resonance_index) + 1;	//note shift
	if (reso_idx_to_be_copied == target_particle_id)
		local_icr2 = 0;
	else
		local_icr2 = lookup_resonance_idx_from_particle_id(reso_idx_to_be_copied) + 1;	//note shift

	try
    {
		Exception::dontPrint();
	
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			hsize_t offset[RANK2D] = {local_icr2 * n_pY_pts + ipY, 0};
			hsize_t count[RANK2D] = {1, chunk_size};
			resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

			// load resonance to be copied first
			resonance_dataset->read(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);

			// now set this to current resonance
			offset[0] = local_icr1 * n_pY_pts + ipY;
			resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
			resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
		}
   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] resonance_chunk;

	return (0);
}

/////////////////////////////////////////////

//int CorrelationFunction::Dump_resonance_HDF_array_spectra(string output_filename, double ******* resonance_array_to_use)
int CorrelationFunction::Dump_resonance_HDF_array_spectra(string output_filename, double * resonance_array_to_use)
{
const int chunk_size = n_pT_pts * n_pphi_pts * /*n_pY_pts **/ qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

	double * resonance_chunk = new double [chunk_size];

	ostringstream filename_stream;
	filename_stream << path << "/" << output_filename;
	ofstream out;
	out.open(filename_stream.str().c_str());

	try
    {
		Exception::dontPrint();
	
		// use loaded chunk to fill resonance_array_to_fill
		for(int ir = 0; ir < NchosenParticle + 1; ir++)
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			hsize_t offset[RANK2D] = {ir * n_pY_pts + ipY, 0};
			hsize_t count[RANK2D] = {1, chunk_size};				// == chunk_dims
			resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
	
			resonance_dataset->read(resonance_array_to_use, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
		}
   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] resonance_chunk;
	out.close();

	return (0);

}

//End of file
