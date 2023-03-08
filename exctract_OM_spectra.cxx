#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>

#include <snfee/snfee.h>
#include <snfee/io/multifile_data_reader.h>

#include <sncabling/om_id.h>
#include <sncabling/gg_cell_id.h>
#include <sncabling/label.h>

#include <snfee/data/raw_event_data.h>
#include <snfee/data/calo_digitized_hit.h>
#include <snfee/data/tracker_digitized_hit.h>

#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMarker.h>

#include <string>

#include "TKEvent/TKEvent/include/TKEvent.h" 

using namespace std;

// these positions of calibration sources migth not perfectly correspond with reality - more precise positions are yet to be determined
double source_ys[6] = {-2087.5, -1252.5, -417.5, 417.5, 1252.5, 2087.5};
double source_zs[7] = {-1425.0, -950.0, -475.0, 0.0, 475.0, 950.0, 1425.0};
double tc_radius    = 22.0;

bool close_to_cal_source(double y, double z, double dist_cut);
double distance_to_OM_front_center(TKOMhit* OM_hit, double x, double y, double z);

int main(int argc, char *argv[])
{
	string input_filename;
	string output_folder;
  
	int run_number = -1;
	int N_events   = -1;

	for (int iarg=1; iarg<argc; ++iarg)
	{
		string arg (argv[iarg]);
	if (arg[0] == '-')
	{
		if (arg=="-i" || arg=="--input")
			input_filename = string(argv[++iarg]);

		else if (arg=="-h" || arg=="--help")
		{
			cout << endl;
			cout << "Usage:   " << argv[0] << " [options]" << endl;
			cout << endl;
			cout << "Options:   -h / --help" << endl;
			cout << "           -i / --input    RED_FILE" << endl;
			cout << "           -r / --run      RUN_NUMBER" << endl;
			cout << "           -o / --output   OUTPUT_FOLDER" << endl;
			cout << "           -n / --n_events NUMBER_OF_EVENTS_TO_READ" << endl;
			cout << endl;
			return 0;
		}
		else if(arg == "-r" || arg == "--run")
		{
			run_number = stod(argv[++iarg]);
		}
		else if(arg == "-o" || arg == "--output")
		{
			output_folder = string(argv[++iarg]);
		}
		else if(arg == "-n" || arg == "--n_events")
		{
			N_events = stoi(argv[++iarg]);
		}
		else
			cerr << "*** unkown option " << arg << endl;
		}
	}

	if (input_filename.empty())
	{
	    cerr << "*** missing input filename !" << endl;
	    return 1;
	}
	if(run_number == -1)
	{
	    cerr << "*** missing run number !" << endl;
	    return 1;
	}
	if(output_folder.empty())
	{
		cerr << "*** missing output folder !" << endl;
		return 1;
	}
	struct stat sb;
  
	if (stat(output_folder.c_str(), &sb) != 0)
	{
		cerr << "*** the path to output folder is invalid !" << endl;
		return 1;
	}

	snfee::initialize();

	/// Configuration for raw data reader
	snfee::io::multifile_data_reader::config_type reader_cfg;
	reader_cfg.filenames.push_back(input_filename);

	// Instantiate a reader
	snfee::io::multifile_data_reader red_source (reader_cfg);

	// Working RED object
	snfee::data::raw_event_data red;
    
	// RED counter
	int red_counter = 0;
  
  	TKEvent *event = new TKEvent();
  
	// cut parameters haven't been optimized yet
	double SOURCE_DIST_CUT = 100.0; // maximum distance between initial electron vertex and calibration source
	int    CHARGE_CUT      = 2000.0; // minimal accepted charge
	double DIST_TRESHOLD   = 200.0; // maximum distance between second electron vertex and center of OM front face
  
	// parameters of charge histograms
	int    N_BINS  = 200;
	double E_MIN   = 0.0;
	double E_MAX   = 3000.0;
	double K_MIN   = 0.035;
	double K_MAX   = 0.075;
	int    K_STEPS = 100;
  
	int N_passed = 0; // number of events accepted by the cut
 
	map<string, vector<TH1F*>> charge_spectra; //OM spectra histograms

	while (red_source.has_record_tag() && (N_events == -1 || red_counter < N_events))
	{
		bool passed = false;
		if(red_counter % 10000 == 0)
			cout << "Event n. " << red_counter << endl;

		// Check the serialization tag of the next record:
		DT_THROW_IF(!red_source.record_tag_is(snfee::data::raw_event_data::SERIAL_TAG),
		          logic_error, "Unexpected record tag '" << red_source.get_record_tag() << "'!");

		// Load the next RED object:
		red_source.load(red);     
	      
		// Load tracker and calo hits so that tracking algorithm can use them
		event = new TKEvent(run_number, red_counter);
		vector<double> OM_charges = {};
		
		const vector<snfee::data::calo_digitized_hit> red_calo_hits = red.get_calo_hits();
		for (const snfee::data::calo_digitized_hit & red_calo_hit : red_calo_hits)
		{
			bool is_HT;
			if(red_calo_hit.is_high_threshold())
			{
				is_HT = true;
			}
			else if(red_calo_hit.is_low_threshold_only())
			{
				is_HT = false;
			}
			else continue;
		
			const snfee::data::timestamp & reference_time = red_calo_hit.get_reference_time();
				
			const sncabling::om_id om_id = red_calo_hit.get_om_id();
			int om_side, om_wall, om_column, om_row, om_num;
			if(om_id.is_main())
			{
				om_side   = om_id.get_side();
				om_column = om_id.get_column();
				om_row    = om_id.get_row();
				om_num = om_side*20*13 + om_column*13 + om_row;
			}

			else if(om_id.is_xwall())
			{
				om_side   = om_id.get_side();
				om_wall   = om_id.get_wall();
				om_column = om_id.get_column();
				om_row    = om_id.get_row();
				om_num = 520 + om_side*64 +  om_wall*32 + om_column*16 + om_row;
			}

			else if(om_id.is_gveto())
			{
				om_side = om_id.get_side();
				om_wall = om_id.get_wall();
				om_column = om_id.get_column();
				om_num = 520 + 128 + om_side*32 + om_wall*16 + om_column;
			}
		
			// Calo-hit is added here
			event->add_OM_hit(om_num, is_HT, reference_time.get_ticks(), red_calo_hit.get_fwmeas_peak_cell());
			OM_charges.push_back(-red_calo_hit.get_fwmeas_charge());
		}

		int side0_gg_hits = 0;
		int side1_gg_hits = 0;
		
		const vector<snfee::data::tracker_digitized_hit> red_tracker_hits = red.get_tracker_hits();
		for (const snfee::data::tracker_digitized_hit & red_tracker_hit : red_tracker_hits)
		{
		  	if(red_tracker_hit.get_cell_id().get_side() == 0)
		  		side0_gg_hits++;
			if(red_tracker_hit.get_cell_id().get_side() == 1)
		  		side1_gg_hits++;
		  	
		  	const sncabling::gg_cell_id gg_id = red_tracker_hit.get_cell_id();
			int srl[3] = {gg_id.get_side(), gg_id.get_row(), gg_id.get_layer()};
			for(const snfee::data::tracker_digitized_hit::gg_times & gg_timestamps : red_tracker_hit.get_times())
			{
				
				int64_t tsp[7];
		
				for (int it = 0; it < 5; it++)
				{
					if(gg_timestamps.get_anode_time(it).get_ticks() != snfee::data::INVALID_TICKS)
					{
						tsp[it] = gg_timestamps.get_anode_time(it).get_ticks();
					}
					else
					{
						tsp[it] = -1;
					}
				}

				if(gg_timestamps.get_bottom_cathode_time().get_ticks() != snfee::data::INVALID_TICKS)
				{
					tsp[5] = gg_timestamps.get_bottom_cathode_time().get_ticks();
				}
				else
				{
					tsp[5] = -1;
				}

				if(gg_timestamps.get_top_cathode_time().get_ticks() != snfee::data::INVALID_TICKS)
				{
					tsp[6] = gg_timestamps.get_top_cathode_time().get_ticks();
				}
				else
				{
					tsp[6] = -1;
				}
		
				// Tracker hit is added here
				event->add_tracker_hit(srl, tsp);
			}
		}

		// calculate tracker hit radii
		event->set_r("Manchester", "distance");
		// reconstructs tracks based on tracker hits and calo-hits (one track for each side)
		event->reconstruct_track(false);
		for(int i = 0;i < event->get_no_tracks();i++)
		{
			TKtrack* track = event->get_track(i);
			
			int track_side = track->get_side();
			
			// check if there are more than 5 hits on given tracker side
			if((track_side == 0 && side0_gg_hits > 5) || (track_side == 1 && side1_gg_hits > 5))
			{
				double a = track->get_a(); // track slope in horizontal plane
				double b = track->get_b(); // track shift in horizontal plane
				double c = track->get_c(); // track slope in vertical plane
				double d = track->get_d(); // track shift in vertical plane
				
				// exclude tracks without z-coordinate
				if(c == 0.0 && d == 0.0)
					break;
				
				// calculate track-OM intersection
				double x,y,z;
			
				x = 435.0;
				if(track_side == 0) 
					x = -x;		
				y = a*x + b;
				
				if(y > 2505.5)
					x = (2505.5 - b) / a;
				else if(y < -2505.5)
					x = (-2505.5 - b) / a;
					
				z = c * x + d;
				if(z > 1550.0)
				{
					x = (1550.0 - d) / c;
					y = a * x + b;
				}
				else if(z < -1550.0)
				{
					x = (-1550.0 - d) / c;
					y = a * x + b;
				}
				
				// check if the first electron vertex is close to a calibration source
				double source_y = -1.0;
				double source_z = -1.0;
				for(int j = 0;j < 6;j++)
				{
					for(int k = 0;k < 7;k++)
					{
						if((source_ys[j] - b) * (source_ys[j] - b) + (source_zs[k] - d) * (source_zs[k] - d) < SOURCE_DIST_CUT*SOURCE_DIST_CUT)
						{
							source_y = source_ys[j];
							source_z = source_zs[k];
						}
					}
				}
				
				if(source_y != -1.0)
				{
					// check if there are tracker hits near the calibration source
					int N_tr_hits = event->get_tr_hits().size();
					int N_close_to_source = 0;
					for(int j = 0;j < N_tr_hits;j++)
					{
						double dist_X = abs(event->get_tr_hit(j)->get_xy('x'));
						double dist_Y = abs(source_y - event->get_tr_hit(j)->get_xy('y'));
						// check if the hit is close to the source
						if(event->get_tr_hit(j)->get_SRL('s') == track_side && dist_X < 11.0 * tc_radius && dist_Y < 10.0 * tc_radius)
							N_close_to_source++;
					}
					if(N_close_to_source < 4)
						break;
					// find the triggered OM closest to the second electron vertex (with measured charge higher than CHARGE_CUT)
					int N_OM_hits = event->get_OM_hits().size();
					double min_distance = 9999999.0;
					double min_charge;
					TKOMhit* min_OM_hit;

		      			for(int j = 0;j < N_OM_hits;j++)
		      			{

			  			double distance = distance_to_OM_front_center(event->get_OM_hit(j), x, y, z);
			  			double charge   = OM_charges.at(j);
						if(distance < min_distance && charge > CHARGE_CUT)
						{
							min_charge   = charge;
							min_distance = distance;
							min_OM_hit   = event->get_OM_hit(j);
						}
					}

					// check if distance between second electron vertex and closest triggered OM is small enough
					if(min_distance < DIST_TRESHOLD)
					{
						// check how many gg cells are in time window of the OM hit
						double OM_time = min_OM_hit->get_OM_TDC() * 6.25E-3; // OM trigger time in us
						int N_hits_in_time_window = 0;
						for(int j = 0;j < N_tr_hits;j++)
						{
							double time_difference = event->get_tr_hit(j)->get_tsp('0') * 12.5E-3 - OM_time;
							if(time_difference > -0.2 && time_difference < 5.0)
								N_hits_in_time_window++;
						}
						if(N_hits_in_time_window > 15)
							break;
							
						int SWCR[4]; // side, wall, collumn, row
						int om_num = min_OM_hit->get_OM_num();
						SWCR[0]    = min_OM_hit->get_SWCR('s');
						SWCR[1]    = min_OM_hit->get_SWCR('w');
						SWCR[2]    = min_OM_hit->get_SWCR('c');
						SWCR[3]    = min_OM_hit->get_SWCR('r');
						
						string om_title;
						
						if(om_num < 520) // MW
							om_title = "1302:" + to_string(SWCR[0]) + "." + to_string(SWCR[2]) + "." + to_string(SWCR[3]) + "_" + to_string(om_num);
						else if(om_num < 684) // xcalo
							om_title = "1232:" + to_string(SWCR[0]) + "." + to_string(SWCR[1]) + "." + to_string(SWCR[2]) + "." + to_string(SWCR[3]) + "_" + to_string(om_num);
						else // gveto
							om_title = "1252:" + to_string(SWCR[0]) + "." + to_string(SWCR[1]) + "." + to_string(SWCR[2]) + "_" + to_string(om_num);
						
						// save measured energy
						if(charge_spectra.find(om_title) == charge_spectra.end())
						{
							charge_spectra[om_title] = {};
							for(int j = 0;j < K_STEPS;j++)
							{
								double k = K_MIN + j * (K_MAX - K_MIN) / K_STEPS;
								charge_spectra[om_title].push_back(new TH1F((om_title + "_k=" + to_string(k)).c_str(), (om_title + "_k=" + to_string(k)).c_str(), N_BINS, E_MIN, E_MAX));
							}
						}
						for(int j = 0;j < K_STEPS;j++)
						{
							double k = K_MIN + j * (K_MAX - K_MIN) / K_STEPS;
							double D = sqrt(x*x + (y - b)*(y - b) + (z - d)*(z - d)); // track length
							// more precise correction will be used in the future
							double E = k * min_charge + (0.251102 * D - 89.787) * exp(-0.00349 * k * min_charge) + 0.0635012 * D - 1.65222;
							charge_spectra[om_title].at(j)->Fill(E);
						}
						passed = true;
					}
				}
			}
		}
		if(passed)
			N_passed++;
		red_counter++;
	}
	
	// save OM spectra
	for(auto item : charge_spectra)
	{
		TFile *file = new TFile((output_folder + item.first + ".root").c_str(), "RECREATE");
		for(int i = 0;i < K_STEPS;i++)
		{
			item.second.at(i)->Write();
		}
  		delete file;
	}
  
	snfee::terminate();
	cout << "Total number of events: " << red_counter << endl;
	cout << "Total number of passed events: " << N_passed << endl;
	cout << "Fraction of passed events: " << (double) N_passed / (double) red_counter << endl;

	return 0;
}

// checks if distance between point with coordinates y, z on the source foil and its closest calibration source is smaler than dist_cut
bool close_to_cal_source(double y, double z, double dist_cut)
{
	for(int i = 0;i < 6;i++)
	{
		for(int j = 0;j < 7;j++)
		{
			if((source_ys[i] - y) * (source_ys[i] - y) + (source_zs[j] - z) * (source_zs[j] - z) < dist_cut*dist_cut)
				return true;
		}
	}
	return false;
}

// calculates the distance between the point with coordinates x, y, z and the front face center of given OM
double distance_to_OM_front_center(TKOMhit* OM_hit, double x, double y, double z)
{
	double OM_coords[3];
	OM_coords[0] = OM_hit->get_xyz('x');
	OM_coords[1] = OM_hit->get_xyz('y');
	OM_coords[2] = OM_hit->get_xyz('z');
	
	int om_num = OM_hit->get_OM_num();
	
	if(om_num < 520) // MW
		return sqrt((y - OM_coords[1]) * (y - OM_coords[1]) + (z - OM_coords[2]) * (z - OM_coords[2]));
	else if(om_num < 684) // xcalo
		return sqrt((x - OM_coords[0]) * (x - OM_coords[0]) + (z - OM_coords[2]) * (z - OM_coords[2]));
	else // gveto
		return sqrt((y - OM_coords[1]) * (y - OM_coords[1]) + (z - OM_coords[2]) * (z - OM_coords[2]));
}
