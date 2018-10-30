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
#include <omp.h>

#include "H5Cpp.h"

#include "cfwr.h"
#include "cfwr_lib.h"
#include "lib.h"

using namespace std;

const int RANK2D = 2;
int return_flag = 0;

void CorrelationFunction::Initialize_HDF_resonance_array()
{
	*global_out_stream_ptr << "Initializing HDF resonance array...";
	int HDFInitializationSuccess = Administrate_resonance_HDF_array(0);
	*global_out_stream_ptr << "...done.\n";

	return;
}

void CorrelationFunction::Close_HDF_resonance_array()
{
	if (!thermal_pions_only)
	{
		*global_out_stream_ptr << "Closing HDF resonance array...";
		int HDFInitializationSuccess = Administrate_resonance_HDF_array(2);
		*global_out_stream_ptr << "...done.\n";
	}

	return;
}

void CorrelationFunction::Reset_HDF_resonance_array()
{
	*global_out_stream_ptr << "Resetting HDF resonance array...";
	int HDFInitializationSuccess = Administrate_resonance_HDF_array(3);
	*global_out_stream_ptr << "...done.\n";

	return;
}

void CorrelationFunction::Initialize_HDF_target_thermal_array()
{
	*global_out_stream_ptr << "Initializing HDF target thermal array...";
	int HDFInitializationSuccess = Administrate_target_thermal_HDF_array(0);
	*global_out_stream_ptr << "...done.\n";

	return;
}

void CorrelationFunction::Initialize_HDF_target_full_array()
{
	*global_out_stream_ptr << "Initializing HDF target full array...";
	int HDFInitializationSuccess = Administrate_target_full_HDF_array(0);
	*global_out_stream_ptr << "...done.\n";

	return;
}

//*******************************************
// HDF array for resonances
//*******************************************
/////////////////////////////////////////////
int CorrelationFunction::Administrate_resonance_HDF_array(int administration_mode)
{
	// administration_mode:
	//	0 - Initialize
	//	1 - Open (already initialized)
	//	2 - Close
	const int chunk_size = n_pT_pts * n_pphi_pts * n_pY_pts * qxnpts * qynpts * ntrig;

	double * resonance_chunk = new double [chunk_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << path << "/resonance_spectra.h5";
	H5std_string RESONANCE_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string RESONANCE_DATASET_NAME("ra");

	try
	{
		//Exception::dontPrint();

		DSetCreatPropList cparms;
		hsize_t dimsm[RANK2D] = {1, chunk_size};
		//hsize_t dims[RANK2D] = { (NchosenParticle + 1) * ((qtnpts+1)/2) * ((qznpts+1)/2), chunk_size};	//only hold quarter of q-space to keep *h5 files as small as possible
		hsize_t dims[RANK2D] = { NchosenParticle + 1, chunk_size};

		switch (administration_mode)
		{
			case 0:	//Initialization
			{
				//bool file_does_not_already_exist = !fexists(filename_stream_ra.str().c_str());
				bool file_does_not_already_exist = true;	//force full initialization for timebeing...

				resonance_file = new H5::H5File(RESONANCE_FILE_NAME, H5F_ACC_TRUNC);

				hsize_t chunk_dims[RANK2D] = {1, chunk_size};
				cparms.setChunk( RANK2D, chunk_dims );

				resonance_dataspace = new H5::DataSpace (RANK2D, dims);

				resonance_dataset = new H5::DataSet( resonance_file->createDataSet(RESONANCE_DATASET_NAME, PredType::NATIVE_DOUBLE, *resonance_dataspace, cparms) );

				hsize_t count[RANK2D] = {1, chunk_size};
				hsize_t offset[RANK2D] = {0, 0};

				resonance_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);

				if (file_does_not_already_exist)
				{
					*global_out_stream_ptr << "HDF resonance file doesn't exist!  Initializing to zero..." << endl;

					for (int ir = 0; ir < NchosenParticle + 1; ++ir)
					//for (int iqt = 0; iqt < (qtnpts+1)/2; ++iqt)
					//for (int iqz = 0; iqz < (qznpts+1)/2; ++iqz)
					{
						//offset[0] = HDF_indexer(ir, iqt, iqz);
						offset[0] = ir;
						resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

						for (int iidx = 0; iidx < chunk_size; ++iidx)
							resonance_chunk[iidx] = 0.0;

						//initialize everything with zeros
						resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
					}
				}
				break;
			}
			case 1:
			{
				resonance_dataspace = new H5::DataSpace (RANK2D, dims);
				resonance_file = new H5::H5File(RESONANCE_FILE_NAME, H5F_ACC_RDWR);
				resonance_dataset = new H5::DataSet( resonance_file->openDataSet( RESONANCE_DATASET_NAME ) );
				resonance_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
				break;
			}
			case 2:
			{
				resonance_memspace->close();
				resonance_dataset->close();
				resonance_file->close();
				delete resonance_memspace;
				delete resonance_file;
				delete resonance_dataset;
				delete resonance_dataspace;
				break;
			}
			case 3:	//new mode for open and reset to zero
			{
				resonance_dataspace = new H5::DataSpace (RANK2D, dims);
				resonance_file = new H5::H5File(RESONANCE_FILE_NAME, H5F_ACC_RDWR);
				resonance_dataset = new H5::DataSet( resonance_file->openDataSet( RESONANCE_DATASET_NAME ) );
				resonance_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
				*global_out_stream_ptr << "Resetting HDF resonance array to zero..." << endl;

				hsize_t count[RANK2D] = {1, chunk_size};
				hsize_t offset[RANK2D] = {0, 0};

				for (int ir = 0; ir < NchosenParticle + 1; ++ir)
				{
					offset[0] = ir;
					resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

					for (int iidx = 0; iidx < chunk_size; ++iidx)
						resonance_chunk[iidx] = 0.0;

					//reset everything to zero
					resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
				}
				break;
			}
			default:
			{
				cerr << __FUNCTION__ << "(): " << "administration_mode = " << administration_mode << " not supported!  Exiting..." << endl;
				exit(1);
				break;
			}
		}
	}

	catch(FileIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "FileIException error!" << endl;
		return_flag = -1;
	}

	catch(H5::DataSetIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSetIException error!" << endl;
		return_flag = -2;
	}

	catch(H5::DataSpaceIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSpaceIException error!" << endl;
		return_flag = -3;
	}

	delete [] resonance_chunk;

	return (return_flag);
}

/////////////////////////////////////////////



int CorrelationFunction::Administrate_target_thermal_HDF_array(int administration_mode)
{
	// administration_mode:
	//	0 - Initialize
	//	1 - Open (already initialized)
	//	2 - Close
	const int chunk_size = n_pT_pts * n_pphi_pts * n_pY_pts * qxnpts * qynpts * ntrig;

	double * tta_chunk = new double [chunk_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << path << "/target_thermal_moments.h5";
	H5std_string TARGET_THERMAL_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string TARGET_THERMAL_DATASET_NAME("tta");

	try
	{
		Exception::dontPrint();

		DSetCreatPropList cparms;
		hsize_t dimsm[RANK2D] = {1, chunk_size};
		hsize_t dims[RANK2D] = {((qtnpts+1)/2) * ((qznpts+1)/2), chunk_size};

		switch (administration_mode)
		{
			case 0:	//Initialization
			{
				//bool file_does_not_already_exist = !fexists(filename_stream_ra.str().c_str());
				bool file_does_not_already_exist = true;	//force full initialization for timebeing...
				tta_file = new H5::H5File(TARGET_THERMAL_FILE_NAME, H5F_ACC_TRUNC);

				//DSetCreatPropList cparms;
				hsize_t chunk_dims[RANK2D] = {1, chunk_size};
				cparms.setChunk( RANK2D, chunk_dims );

				//hsize_t dims[RANK2D] = {((qtnpts+1)/2) * ((qznpts+1)/2), chunk_size};
				tta_dataspace = new H5::DataSpace (RANK2D, dims);

				tta_dataset = new H5::DataSet( tta_file->createDataSet(TARGET_THERMAL_DATASET_NAME, PredType::NATIVE_DOUBLE, *tta_dataspace, cparms) );

				hsize_t count[RANK2D] = {1, chunk_size};
				//hsize_t dimsm[RANK2D] = {1, chunk_size};
				hsize_t offset[RANK2D] = {0, 0};

				tta_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
				if (file_does_not_already_exist)
				{
					*global_out_stream_ptr << "HDF thermal target moments file doesn't exist!  Initializing to zero..." << endl;

					for (int iqt = 0; iqt < (qtnpts+1)/2; ++iqt)
					for (int iqz = 0; iqz < (qznpts+1)/2; ++iqz)
					{
						offset[0] = HDF_indexer(0, iqt, iqz);
						tta_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

						for (int iidx = 0; iidx < chunk_size; ++iidx)
							tta_chunk[iidx] = 0.0;

						//initialize everything with zeros
						tta_dataset->write(tta_chunk, PredType::NATIVE_DOUBLE, *tta_memspace, *tta_dataspace);
					}
				}
				break;
			}
			case 1:
			{
				//hsize_t dimsm[RANK2D] = {1, chunk_size};
				//hsize_t dims[RANK2D] = {((qtnpts+1)/2) * ((qznpts+1)/2), chunk_size};

				tta_dataspace = new H5::DataSpace (RANK2D, dims);
				tta_file = new H5::H5File(TARGET_THERMAL_FILE_NAME, H5F_ACC_RDWR);
				tta_dataset = new H5::DataSet( tta_file->openDataSet( TARGET_THERMAL_DATASET_NAME ) );
				tta_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
				break;
			}
			case 2:
			{
				tta_memspace->close();
				tta_dataset->close();
				tta_file->close();
				delete tta_memspace;
				delete tta_file;
				delete tta_dataset;
				delete tta_dataspace;
				break;
			}
			default:
			{
				cerr << __FUNCTION__ << "(): " << "administration_mode = " << administration_mode << " not supported!  Exiting..." << endl;
				exit(1);
				break;
			}
		}
	}

	catch(FileIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "FileIException error!" << endl;
		return_flag = -1;
	}

	catch(H5::DataSetIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSetIException error!" << endl;
		return_flag = -2;
	}

	catch(H5::DataSpaceIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSpaceIException error!" << endl;
		return_flag = -3;
	}

	delete [] tta_chunk;

	return (return_flag);
}

/////////////////////////////////////////////

int CorrelationFunction::Administrate_target_full_HDF_array(int administration_mode)
{
	// administration_mode:
	//	0 - Initialize
	//	1 - Open (already initialized)
	//	2 - Close
	const int chunk_size = n_pT_pts * n_pphi_pts * n_pY_pts * qxnpts * qynpts * ntrig;

	double * tfa_chunk = new double [chunk_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << path << "/target_full_moments.h5";
	H5std_string TARGET_FULL_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string TARGET_FULL_DATASET_NAME("tfa");

	try
	{
		Exception::dontPrint();

		DSetCreatPropList cparms;
		hsize_t dimsm[RANK2D] = {1, chunk_size};
		hsize_t dims[RANK2D] = {((qtnpts+1)/2) * ((qznpts+1)/2), chunk_size};

		switch (administration_mode)
		{
			case 0:	//Initialization
			{
				bool file_does_not_already_exist = true;
				tfa_file = new H5::H5File(TARGET_FULL_FILE_NAME, H5F_ACC_TRUNC);

				hsize_t chunk_dims[RANK2D] = {1, chunk_size};
				cparms.setChunk( RANK2D, chunk_dims );

				tfa_dataspace = new H5::DataSpace (RANK2D, dims);

				tfa_dataset = new H5::DataSet( tfa_file->createDataSet(TARGET_FULL_DATASET_NAME, PredType::NATIVE_DOUBLE, *tfa_dataspace, cparms) );

				hsize_t count[RANK2D] = {1, chunk_size};
				hsize_t offset[RANK2D] = {0, 0};

				tfa_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
				if (file_does_not_already_exist)
				{
					*global_out_stream_ptr << "HDF full target moments file doesn't exist!  Initializing to zero..." << endl;

					for (int iqt = 0; iqt < (qtnpts+1)/2; ++iqt)
					for (int iqz = 0; iqz < (qznpts+1)/2; ++iqz)
					{
						offset[0] = HDF_indexer(0, iqt, iqz);
						tfa_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

						for (int iidx = 0; iidx < chunk_size; ++iidx)
							tfa_chunk[iidx] = 0.0;

						//initialize everything with zeros
						tfa_dataset->write(tfa_chunk, PredType::NATIVE_DOUBLE, *tfa_memspace, *tfa_dataspace);
					}
				}
				break;
			}
			case 1:
			{
				tfa_dataspace = new H5::DataSpace (RANK2D, dims);
				tfa_file = new H5::H5File(TARGET_FULL_FILE_NAME, H5F_ACC_RDWR);
				tfa_dataset = new H5::DataSet( tfa_file->openDataSet( TARGET_FULL_DATASET_NAME ) );
				tfa_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
				break;
			}
			case 2:
			{
				tfa_memspace->close();
				tfa_dataset->close();
				tfa_file->close();
				delete tfa_memspace;
				delete tfa_file;
				delete tfa_dataset;
				delete tfa_dataspace;
				break;
			}
			default:
			{
				cerr << __FUNCTION__ << "(): " << "administration_mode = " << administration_mode << " not supported!  Exiting..." << endl;
				exit(1);
				break;
			}
		}
	}

	catch(FileIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "FileIException error!" << endl;
		return_flag = -1;
	}

	catch(H5::DataSetIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSetIException error!" << endl;
		return_flag = -2;
	}

	catch(H5::DataSpaceIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSpaceIException error!" << endl;
		return_flag = -3;
	}

	delete [] tfa_chunk;

	return (return_flag);
}

/////////////////////////////////////////////

int CorrelationFunction::Access_resonance_in_HDF_array(int local_pid, int iqt, int iqz, int access_mode, double * resonance_array_to_use, bool verbose)
{
	// access_mode:
	//	0 - set array chunk
	//	1 - get array chunk
	const int chunk_size = n_pT_pts * n_pphi_pts * n_pY_pts * qxnpts * qynpts * ntrig;

	int local_icr;
	if (local_pid == target_particle_id)
		local_icr = 0;
	else
		local_icr = lookup_resonance_idx_from_particle_id(local_pid) + 1;	//note shift

	try
	{
		Exception::dontPrint();
		//hsize_t offset[RANK2D] = {HDF_indexer(local_icr, iqt, iqz), 0};
		hsize_t offset[RANK2D] = {local_icr, 0};
		//if (verbose) cout << "In Access_resonance_in_HDF_array(...): called with/using arguments " << local_pid << "   " << iqt << "   " << iqz << "   " << access_mode << "   " << local_icr << "   " << HDF_indexer(local_icr, iqt, iqz) << endl;
		hsize_t count[RANK2D] = {1, chunk_size};				// == chunk_dims
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
		//if (verbose) debugger(__LINE__, __FILE__);
	
		switch(access_mode)
		{
			case 0:
			{
				//if (verbose) debugger(__LINE__, __FILE__);
				resonance_dataset->write(resonance_array_to_use, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
				break;
			}
			case 1:
			{
				//if (verbose) debugger(__LINE__, __FILE__);
				resonance_dataset->read(resonance_array_to_use, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
				break;
			}
			default:
			{
				cerr << __FUNCTION__ << "(): " << "access_mode = " << access_mode << " not supported!  Exiting..." << endl;
				exit(1);
				break;
			}
		}
   }

	catch(FileIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "FileIException error!" << endl;
		return_flag = -1;
	}

	catch(H5::DataSetIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSetIException error!" << endl;
		return_flag = -2;
	}

	catch(H5::DataSpaceIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSpaceIException error!" << endl;
		return_flag = -3;
	}

	//if (verbose) debugger(__LINE__, __FILE__);

	return (return_flag);
}

/////////////////////////////////////////////

int CorrelationFunction::Access_target_thermal_in_HDF_array(int iqt, int iqz, int access_mode, double * tta_array_to_use, bool verbose)
{
	// access_mode:
	//	0 - set array chunk
	//	1 - get array chunk
	const int chunk_size = n_pT_pts * n_pphi_pts * n_pY_pts * qxnpts * qynpts * ntrig;

	try
	{
		Exception::dontPrint();

		hsize_t offset[RANK2D] = {HDF_indexer(0, iqt, iqz), 0};
		hsize_t count[RANK2D] = {1, chunk_size};				// == chunk_dims
		tta_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
		//if (verbose) cout << "In Access_target_thermal_in_HDF_array(...): called with/using arguments " << iqt << "   " << iqz << "   " << access_mode << "   " << 0 << "   "
		//					<< chunk_size << "   " << n_pT_pts * n_pphi_pts * n_pY_pts * qxnpts * qynpts * ntrig << "   " << HDF_indexer(0, iqt, iqz) << endl;

		switch(access_mode)
		{
			case 0:
			{
				//debugger(__LINE__, __FILE__);
				tta_dataset->write(tta_array_to_use, PredType::NATIVE_DOUBLE, *tta_memspace, *tta_dataspace);
				break;
			}
			case 1:
			{
				tta_dataset->read(tta_array_to_use, PredType::NATIVE_DOUBLE, *tta_memspace, *tta_dataspace);
				break;
			}
			default:
			{
				cerr << __FUNCTION__ << "(): " << "access_mode = " << access_mode << " not supported!  Exiting..." << endl;
				exit(1);
				break;
			}
		}

   }

	catch(FileIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "FileIException error!" << endl;
		return_flag = -1;
	}

	catch(H5::DataSetIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSetIException error!" << endl;
		return_flag = -2;
	}

	catch(H5::DataSpaceIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSpaceIException error!" << endl;
		return_flag = -3;
	}

	return (return_flag);
}



/////////////////////////////////////////////

int CorrelationFunction::Access_target_full_in_HDF_array(int iqt, int iqz, int access_mode, double * tfa_array_to_use, bool verbose)
{
	// access_mode:
	//	0 - set array chunk
	//	1 - get array chunk
	const int chunk_size = n_pT_pts * n_pphi_pts * n_pY_pts * qxnpts * qynpts * ntrig;

	try
	{
		Exception::dontPrint();

		hsize_t offset[RANK2D] = {HDF_indexer(0, iqt, iqz), 0};
		hsize_t count[RANK2D] = {1, chunk_size};				// == chunk_dims
		tfa_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
		//if (verbose)
		//	cout << "In Access_full_thermal_in_HDF_array(...): called with/using arguments "
		//			<< iqt << "   " << iqz << "   " << access_mode << "   " << 0 << "   "
		//			<< chunk_size << "   " << n_pT_pts * n_pphi_pts * n_pY_pts * qxnpts * qynpts * ntrig
		//			<< "   " << HDF_indexer(0, iqt, iqz) << endl;

		switch(access_mode)
		{
			case 0:
			{
				tfa_dataset->write(tfa_array_to_use, PredType::NATIVE_DOUBLE, *tfa_memspace, *tfa_dataspace);
				break;
			}
			case 1:
			{
				tfa_dataset->read(tfa_array_to_use, PredType::NATIVE_DOUBLE, *tfa_memspace, *tfa_dataspace);
				break;
			}
			default:
			{
				cerr << __FUNCTION__ << "(): " << "access_mode = " << access_mode << " not supported!  Exiting..." << endl;
				exit(1);
				break;
			}
		}

   }

	catch(FileIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "FileIException error!" << endl;
		return_flag = -1;
	}

	catch(H5::DataSetIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSetIException error!" << endl;
		return_flag = -2;
	}

	catch(H5::DataSpaceIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSpaceIException error!" << endl;
		return_flag = -3;
	}

	return (return_flag);
}

/////////////////////////////////////////////

int CorrelationFunction::Copy_chunk(int current_resonance_index, int reso_idx_to_be_copied)
{
	const int chunk_size = n_pT_pts * n_pphi_pts * n_pY_pts * qxnpts * qynpts * ntrig;

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

		//for (int iqt = 0; iqt < (qtnpts+1)/2; ++iqt)
		//for (int iqz = 0; iqz < (qznpts+1)/2; ++iqz)
		//{
			//hsize_t offset[RANK2D] = {HDF_indexer(local_icr2, iqt, iqz), 0};
			hsize_t offset[RANK2D] = {local_icr2, 0};
			hsize_t count[RANK2D] = {1, chunk_size};
			resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

			// load resonance to be copied first
			resonance_dataset->read(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);

			// now set this to current resonance
			//offset[0] = HDF_indexer(local_icr1, iqt, iqz);
			offset[0] = local_icr1;
			resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
			resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
		//}
   }

	catch(FileIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "FileIException error!" << endl;
		return_flag = -1;
	}

	catch(H5::DataSetIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSetIException error!" << endl;
		return_flag = -2;
	}

	catch(H5::DataSpaceIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSpaceIException error!" << endl;
		return_flag = -3;
	}

	delete [] resonance_chunk;

	return (return_flag);
}


//////////////////////////////////////////////////////////
// Functions to read and write Bessel coefficient arrays

int CorrelationFunction::Administrate_besselcoeffs_HDF_array(int administration_mode, int particle_mode /*==0*/)
{
	// administration_mode:
	//	0 - Initialize
	//	1 - Open (already initialized)
	//	2 - Close
	// particle_mode:
	//	0 - Heavy resonances (default)
	//	1 - Light particles (pions) --> need different range and n_alpha_points
	const int n_chunks = n_pY_pts;
	//const int chunk_size = 4 * FO_length * n_alpha_points;
	int chunk_size = 4 * FO_length * n_alpha_points;
	//if (particle_mode)
	//	chunk_size = 4 * FO_length * n_alpha_points_PIONS;

	double * besselcoeffs_chunk = new double [chunk_size];

	string pions_stem = "";
	//if (particle_mode)
	//	pions_stem = "_PIONS";

	ostringstream filename_stream_ra;
	filename_stream_ra << path << "/Bessel_coefficients" << pions_stem << ".h5";
	H5std_string BC_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string BC_DATASET_NAME("bc");

	try
	{
		Exception::dontPrint();

		DSetCreatPropList cparms;
		hsize_t dimsm[RANK2D] = {1, chunk_size};
		hsize_t dims[RANK2D] = { n_chunks, chunk_size};

		switch (administration_mode)
		{
			case 0:	//Initialization
			{
				bool file_does_not_already_exist = !fexists(filename_stream_ra.str().c_str());
				//bool file_does_not_already_exist = true;	//force full initialization for timebeing...
				besselcoeffs_file = new H5::H5File(BC_FILE_NAME, H5F_ACC_TRUNC);

				//DSetCreatPropList cparms;
				hsize_t chunk_dims[RANK2D] = {1, chunk_size};
				cparms.setChunk( RANK2D, chunk_dims );

				//hsize_t dims[RANK2D] = {n_chunks, chunk_size};
				besselcoeffs_dataspace = new H5::DataSpace (RANK2D, dims);

				besselcoeffs_dataset = new H5::DataSet( besselcoeffs_file->createDataSet(BC_DATASET_NAME, PredType::NATIVE_DOUBLE, *besselcoeffs_dataspace, cparms) );

				hsize_t count[RANK2D] = {1, chunk_size};
				//hsize_t dimsm[RANK2D] = {1, chunk_size};
				hsize_t offset[RANK2D] = {0, 0};

				besselcoeffs_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
				if (file_does_not_already_exist)
				{
					*global_out_stream_ptr << "HDF besselcoeffs file doesn't exist!  Initializing to zero..." << endl;

					for (int ipY = 0; ipY < n_pY_pts; ++ipY)
					{
						offset[0] = ipY;
						besselcoeffs_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

						for (int iidx = 0; iidx < chunk_size; ++iidx)
							besselcoeffs_chunk[iidx] = 0.0;

						//initialize everything with zeros
						besselcoeffs_dataset->write(besselcoeffs_chunk, PredType::NATIVE_DOUBLE, *besselcoeffs_memspace, *besselcoeffs_dataspace);
					}
				}
				break;
			}
			case 1:	//Open
			{
				//hsize_t dimsm[RANK2D] = {1, chunk_size};
				//hsize_t dims[RANK2D] = { n_chunks, chunk_size};

				besselcoeffs_dataspace = new H5::DataSpace (RANK2D, dims);
				besselcoeffs_file = new H5::H5File(BC_FILE_NAME, H5F_ACC_RDWR);
				besselcoeffs_dataset = new H5::DataSet( besselcoeffs_file->openDataSet( BC_DATASET_NAME ) );
				besselcoeffs_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
				break;
			}
			case 2:	//Close
			{
				besselcoeffs_memspace->close();
				besselcoeffs_dataset->close();
				besselcoeffs_file->close();
				delete besselcoeffs_memspace;
				delete besselcoeffs_file;
				delete besselcoeffs_dataset;
				delete besselcoeffs_dataspace;
				break;
			}
			default:
			{
				cerr << __FUNCTION__ << "(): " << "administration_mode = " << administration_mode << " not supported!  Exiting..." << endl;
				exit(1);
				break;
			}
		}
	}

	catch(FileIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "FileIException error!" << endl;
		return_flag = -1;
	}

	catch(H5::DataSetIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSetIException error!" << endl;
		return_flag = -2;
	}

	catch(H5::DataSpaceIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSpaceIException error!" << endl;
		return_flag = -3;
	}

	delete [] besselcoeffs_chunk;

	return (return_flag);
}

/////////////////////////////////////////////

int CorrelationFunction::Access_besselcoeffs_in_HDF_array(int ipY, int access_mode, double * besselcoeffs_array_to_use, int particle_mode /*==0*/)
{
	// access_mode:
	//	0 - set array chunk
	//	1 - get array chunk
	// particle_mode:
	//	0 - Heavy resonances (default)
	//	1 - Light particles (pions) --> need different range and n_alpha_points
	const int n_chunks = n_pY_pts;
	//const int chunk_size = 4 * FO_length * n_alpha_points;
	int chunk_size = 4 * FO_length * n_alpha_points;
	//if (particle_mode)
	//	chunk_size = 4 * FO_length * n_alpha_points_PIONS;

	double * besselcoeffs_chunk = new double [chunk_size];

	try
	{
		Exception::dontPrint();
		hsize_t offset[RANK2D] = {ipY, 0};
		hsize_t count[RANK2D] = {1, chunk_size};				// == chunk_dims
		besselcoeffs_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		switch(access_mode)
		{
			case 0:
			{
				besselcoeffs_dataset->write(besselcoeffs_array_to_use, PredType::NATIVE_DOUBLE, *besselcoeffs_memspace, *besselcoeffs_dataspace);
				break;
			}
			case 1:
			{
				besselcoeffs_dataset->read(besselcoeffs_array_to_use, PredType::NATIVE_DOUBLE, *besselcoeffs_memspace, *besselcoeffs_dataspace);
				break;
			}
			default:
			{
				cerr << __FUNCTION__ << "(): " << "access_mode = " << access_mode << " not supported!  Exiting..." << endl;
				exit(1);
				break;
			}
		}
   }

	catch(FileIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "FileIException error!" << endl;
		return_flag = -1;
	}

	catch(H5::DataSetIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSetIException error!" << endl;
		return_flag = -2;
	}

	catch(H5::DataSpaceIException error)
	{
		error.printError();
		cerr << __FUNCTION__ << "(): " << "DataSpaceIException error!" << endl;
		return_flag = -3;
	}

	delete [] besselcoeffs_chunk;

	return (return_flag);
}

//End of file
