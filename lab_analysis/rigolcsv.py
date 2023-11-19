from pandas import DataFrame, Series
from numpy import *
import csv

class RigolCSV:
    """ Class for clean handling Rigol spectrum analyzer data """
    def __init__(self, filename, exclude=2, trace2=False, no_offset=True, zero=True):
        """
        'exclude': the number of column to start at; i.e. exclude=2 means omit column 1.#TODO rename
        'trace2': True if two traces. False by default. #TODO generalize for N traces
        'no_offset': if True shifts xpts so maximum in y is at x=0. True by default.
        'zero': set the height of each trace to 0dB
        
        ##TODO: need much better way of linking trace with corresponding data.. maybe just append the 
        points lists as a tuple in a list of traces, where trace is no longer a dataframe, which is 
        unnecessary.
        """
        
        self.filename = filename
        self.data = DataFrame(self.process_csv(self.filename, exclude)) # a pandas DataFrame
        self.no_offset = no_offset
        self.zero = zero
        self.xpts2 = None
        self.ypts2 = None
        
        # data from the first trace
        self.trace1_df = self.data[[0,2]].astype(float)
        self.xpts1, self.ypts1 = self.trace1_df.values.transpose()
        self.x0 = 0
        
        # data from second trace
        if trace2 != False:
            self.trace2_df = self.data[[4,6]].astype(float)
            self.xpts2, self.ypts2 = self.trace2_df.values.transpose()
            
        if self.zero == True:
            self.ypts1 -= max(self.ypts1)
            if self.xpts2 is not None:
                self.ypts2 -= max(self.ypts2)
            
        if self.no_offset == True:
            
            if self.xpts2 is not None:
                self.x0 = self.xpts2[where(self.ypts2 == max(self.ypts2))[0]][0]
                self.xpts2 -= self.x0
            else:
                self.x0 = self.xpts1[where(self.ypts1 == max(self.ypts1))[0]][0]
           
            self.xpts1 -= self.x0
        
    def process_csv(self, filename, exclude=2):
        """ exclude value is the row integer index we start with; i.e. exclude=2 leaves out row 1. """
        exampleFile = open(filename, encoding='utf-8')
        exampleReader = csv.reader(exampleFile)
        exampleData = list(exampleReader)
        exampleFile.close()
        return exampleData[exclude:]
    
    def get_region_df(self, offset, halfwidth):
        """ return subset of rows zoomed in on the region within 
            [offset-halfwidth, offset+halfwidth]. 
        """
        sub_dfs = []
        for i, trace in enumerate(self.traces):
            sub_dfs.append(trace[(trace[4*i] > offset - halfwidth) & 
                                 (trace[4*i] < offset + halfwidth)])
        return sub_dfs
        
    @staticmethod
    def get_linewidth(xpts, ypts, ystep=10):
        """for xpts and ypts from a spectral trace with ypts in [dB],
        return linewidth in kHz
        
        ystep can be set to ensure the right-side 3dB point does not
        fall on the left, adjacent to the left-side 3dB point by
        accident. 
        """
        y0 = max(ypts)
        # print(y0)
        idxL = idxR = 0
        for i in range(len(ypts)):
            if ypts[i] > y0 - 3:
                idxL = i
                for j in range(idxL+ystep,len(ypts)):
                    if ypts[j] < ypts[idxL]:
                        idxR = j
                        # print(f"xpts: {xpts[idxL]/1e6,xpts[idxR]/1e6}")
                        # print(f"ypts: {ypts[idxL],ypts[idxR]}")
                        break
                break
        return (xpts[idxR] - xpts[idxL])/1e3
