void createWeightFile(TString inputFile, TString outputFile){
 
 std::ifstream indata;
 indata.open(inputFile.Data());
 if(!indata) {
  std::cerr << "Error: file " << inputFile.Data() << " could not be opened" << std::endl;
  return 1;
 }
 
 ofstream myfile;
 myfile.open (outputFile.Data());
 myfile << "double ptHWeight(double pt){ " << std::endl;
 
 myfile << "if (pt < 0) return 1;" << std::endl;
 
 double x1;
 double x2;
 double y;
 
 
 while ( !indata.eof() ) { // keep reading until end-of-file
  indata >> x1;
  indata >> x2;
  indata >> y;

  myfile << "if (pt >= " << x1 << " && pt < " << x2 << ") return " << y << " ; " << std::endl;  
 } 
 
 myfile << "return " << y << " ; " << std::endl;  
 
 myfile << "} " << std::endl;
 myfile.close(); 
 
 
 ///==== exit ====
 gApplication->Terminate(0);
 
}