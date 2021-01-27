#include <FastQFile.h>
#include <FastQStatus.h>
#include "BaseAsciiMap.h"
#include <stdio.h>

int main(int argc, char ** argv)
{
   // Check for the appropriate number of arguments.
   if(argc != 2)
   {
      printf("./a.out <inputFile>\n");
      exit(-1);
   }

   FastQFile fastQFile(4,4);
   String filename = argv[1];
   // Open the fastqfile with the default UNKNOWN space type which will determine the 
   // base type from the first character in the sequence.
   if(fastQFile.openFile(filename, BaseAsciiMap::UNKNOWN) != FastQStatus::FASTQ_SUCCESS)
   {
      // Failed to open the specified file.
      perror("Failed to open file:");
      return (-1);
   }
   // Keep reading the file until there are no more fastq sequences to process.
   while (fastQFile.keepReadingFile())
   {
      // Read one sequence. This call will read all the lines for 
      // one sequence.
      /////////////////////////////////////////////////////////////////
      // NOTE: It is up to you if you want to process only for success:
      //    if(readFastQSequence() == FASTQ_SUCCESS)
      // or for FASTQ_SUCCESS and FASTQ_INVALID: 
      //    if(readFastQSequence() != FASTQ_FAILURE)
      // Do NOT try to process on a FASTQ_FAILURE
      /////////////////////////////////////////////////////////////////
      if(fastQFile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS)
      {
         // The sequence is valid.
         // For example if you want to print the lines of the sequence:
        // printf("The Sequence ID Line is: %s", fastQFile.mySequenceIdLine.c_str());
         printf("The Sequence ID is: %s\n", fastQFile.mySequenceIdentifier.c_str());
         //printf("The Sequence Line is: %s", fastQFile.myRawSequence.c_str());
         //printf("The Plus Line is: %s", fastQFile.myPlusLine.c_str());
         //printf("The Quality String Line is: %s", fastQFile.myQualityString.c_str());
      }
   }
   // Finished processing all of the sequences in the file.
   // Close the input file.
   fastQFile.closeFile();
   return 0; // It is up to you to determine your return.


 }
