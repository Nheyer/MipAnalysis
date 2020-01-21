#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h> 
#define DEBUG false


std::string get_before_space(std::string full , std::string& left){
    std::string rt = "";
    bool cleaning = false;
#if DEBUG
    std::cerr <<  full << "\t" << full.length() << std::endl;
#endif 
    for(int i = 0; i < full.length(); ++i){
#if DEBUG
        std::cerr << full[i] << " " << ( full[i] == ' ' || full[i] == '\t');
#endif
        if( full[i] == ' ' || full[i] == '\t'){
            cleaning = true;
        } else if (not cleaning){
            rt += full[i];
        } else if (cleaning){
            left += full[i];
        }
    }
    return rt;
}


class read {
    public:
        std::string seq;
        std::string name;
        int qr_tag = -1 ;
        std::string bx_tag = "" ;
        std::string bx_qual = "";
        std::string q_arr;
        std::string remaining = "";
        read(std::string in_seq, std::string in_name, std::string  in_qs){
            //std::cerr << "In functiin to make read " <<  std::endl;
            seq = in_seq;
            //std::cerr << "Assigned the sequence " << std::endl;
            name = get_before_space(in_name, remaining);
            //std::cerr << "Assigned name, now for qual arr " << std::endl;
            q_arr = in_qs;
        };

        void sum_line(){
            long unsigned int sum_q = 0 ;
            //std::cerr << "Sequence \"" << seq << "\" length is " << seq.length() << std::endl;
            for (int i = 0; i < seq.length(); ++i) {
                //std::cerr << "bp " << i <<  " has a quality of \n";// << q_arr[i]  << std::endl;
                sum_q += (int) q_arr[i] - 33 ;  // illumina offset
            }
             qr_tag = sum_q / seq.length();
        }
        void get_bx(int x , int y,read& r2){
        if(x < 0){
            return;
        } else {
            // check they have same name till a "/" or a " "
            for (int i = 0; i < name.length() ; ++i) {
                if(name[i] != r2.name[i]){
                    std::cerr << "Names: " << name << " and " << r2.name << " do not match!\n";
                    exit(-1);
                } else if( name[i] == ' ' or name[i] == '/'){
                    break;
                }
            }
        }
        bx_tag = seq.substr(0, x) + "-MIP-" + r2.seq.substr(0, y);
        r2.bx_tag = bx_tag;
        seq = seq.substr(x, seq.length());
        r2.seq = r2.seq.substr(y, r2.seq.length());
        bx_qual = q_arr.substr(0, x) + "-MIP-" + r2.q_arr.substr(0, y);
        r2.bx_qual = bx_qual;
        q_arr = q_arr.substr(y, q_arr.length());
        r2.q_arr = r2.q_arr.substr(y, r2.q_arr.length());
        }
        void get_read(std::ofstream & out){
            out << name;// << "\t" << remaining;
            if (qr_tag >= 0){
                out << "\t" << "RQ:i:" << qr_tag;
            }
            if (bx_tag != "" && bx_tag.length() == bx_qual.length()){
                out << "\t" << "BX:Z:" << bx_tag 
                    << "\t" << "BQ:Z:" << bx_qual ;
            } else if (bx_tag != ""){
                std::cerr << "We had a BX but no quality, will attempt to output, but it may be malformed!!\n";
                out << "\t" << "BX:Z:" << bx_tag         
                    << "\t" << "BQ:Z:" << bx_qual;
            }
            out << "\n" <<  seq + "\n"
            << "+\n"
            << q_arr << std::endl;
        }
};


int main(int argc, char ** argv) {
    if(argc < 3) {
        std::cerr << "Usage: fastRQ input1.fastq input2.fastq [output_modification] <int> <int> \n"
                  << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
                  << "first two inputs should be the fastq files of read one and read two\n"
                  << "then the modification you wish to be added to the file\n"
                  << "and then finally two ints for bx length on R1 and R2 respectively\n";
    }
    std::cout << "Obtained " << argc << " inputs\n";
    std::string name_line[2] = {};
    std::string seq_line[2] = {};
    std::string q_line[2] = {};
    std::string check[2] = {};
    std::ifstream input[2] ;
    std::ofstream output[2];
    bool interleaf = false;
    bool open ;
    std::cout << "Atempthing to open " << argv[1] << " and " <<  argv[2] << " as inputs\n" ;
    input[0].open(argv[1]);
    input[1].open(argv[2]);
    std::cout << "Opened file " << argv[1] << " and " << argv[2] << " for input\n";
    output[0].open(((std::string) argv[1] + (std::string) argv[3]).c_str());
    if(argc > 5){
        if(argv[5][0] == 't'){
            interleaf = true;
        }
    }
    if (not interleaf) {
        output[1].open(((std::string) argv[2] + (std::string) argv[3]).c_str());
        open = output[1].good();
    } else {
        open = true;
    }
    long long unsigned int count = 0;
    int r1_bx_length = -1;
    int r2_bx_length = -1;
    if(argc > 4){
        r1_bx_length = atoi(argv[4]);
        r2_bx_length = atoi(argv[5]);
    }
    if(input[0].good() && input[1].good() && output[0].good() && open ){
        std::cerr << "All Files were opened as expected\n" ;
    } else {
        std::cerr << "couldn't open file\n";
        exit(-1);
    }
    while(!input[0].eof() && !input[1].eof()){
        for (int i = 0; i < 2 ; ++i) {
            check[i] = "";
            std::getline(input[i], name_line[i]);
            std::getline(input[i], seq_line[i]);
            std::getline(input[i], check[i]);
            std::getline(input[i], q_line[i]);
        }
#if DEBUG
        std::cerr << name_line << std::endl
                  << seq_line  << std::endl
                  << check     << std::endl
                  << q_line    << std::endl;
#endif
        if (check[0] == "+" && check[1] == "+"){
#if DEBUG
            std::cerr << "We can read it, making a read object\n";
#endif
            read read1(seq_line[0], name_line[0], q_line[0]);
            read read2(seq_line[1], name_line[1], q_line[1]);
#if DEBUG
            std::cerr << "suming qulaity lines\n";
#endif
            read1.get_bx(r1_bx_length, r2_bx_length, read2);

#if DEBUG
            std::cerr << "riping off " << x << "bases as our barcode!";
#endif
            read1.sum_line();
            read2.sum_line();
#if DEBUG
            std::cerr << "printing output\n";
#endif
            read1.get_read(output[0]);
            if (interleaf){
                read2.get_read(output[0]);
            } else {
                read2.get_read(output[1]);
            }
            count++;
        }
    }
    // close ALL THE FILES
    std::cout << "Obtained <" << count << "> reads... \n Now closing files!\n";
    input[0].close();
    input[1].close();
    output[0].close();
    if (not interleaf) {
        output[1].close();
    }
}
