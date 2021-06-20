import pandas as pd
from datetime import datetime
from matplotlib import pyplot as plt
#####################################
#Transfer to mgf Class
class MGF:
    def __init__(self):
        self.spectra = []
        self.filename = None

        #filter attributes
        self.filter_peaks = None
        self.pos_spectra = []
        self.neg_spectra = []

        ##difference attributes (Tree)
        #self.diff_root = None
        #self.diff_count = 0
        
        #difference attributes (List)
        self.diff_list = []
        self.unique_diffs = []
        self.diff_counts = []
        self.diff_spectra = []


    def add_spectra(self, spectra):
        self.spectra += spectra

    def save(self, filename = False):
        """saves a dictionary generated from a.mgf file back into a .mgf file

        Keyword Arguments:
        mgf -- dictionary generated from a .mgf file
        filename -- filename to give the new .mgf file
        """
        if filename == False:
            filename_mgf = self.filename + ("_PROCESSED.mgf")
        else:
            filename_mgf = str(filename) + ".mgf"
        #Write to new mgf file
        #Blank file
        with open(filename_mgf, "w") as f:
            f.write("")
        new_file_write = open("filename_mgf", "w")
        new_file_write.write("")
        #Append each line to file
        for entry in the_file:
            with open(filename_mgf, "ab") as f:
                f.write(bytearray(b"BEGIN IONS\n"))
                f.write(bytearray("TITLE=index=" + entry.title + "\n", "utf-8"))
                f.write(bytearray("PEPMASS=" + entry.pepmass + "\n", "utf-8"))
                f.write(bytearray("CHARGE=" + entry.charge + "+\n", "utf-8"))
                f.write(bytearray("SCANS=" + entry.scans + "\n", "utf-8"))
                f.write(bytearray("RTINSECONDS=" + entry.rtinseconds + "\n", "utf-8"))
                for peak in entry.spectra:
                    f.write(bytearray(str(peak[0]) + " " + str(peak[1]) + "\n", "utf-8"))
                f.write(b"END IONS\n\n")

    def filter(self, mz_list, intensity_tol = 0.2, wobble_tol = 0.01, checks = 1, negative = False):
        """Takes a .mgf file converted into a dictionary by extract_mgf() and extracts all spectra with peaks of      a desired m/z value and 
        intensity relative to the max intensity in the spectra
        
        Keyword Arguments:
        mgf_list -- dictionary generated from a .mgf file
        mz_list -- a list of desired m/z values to check spectra for
        intensity_tol -- intensity relative to the max peak in each spectra that a peak must be above to be           detected
        wobble_tol -- value above and below the specified m/z values peaks can be detected
        checks -- the number of peaks that a spectra must have with the specified m/z value to be detected, a value below 1 will require 1 or more peaks to be present and a value greater than the amount of specified peaks will require all peaks to be present
        """
        #Extract spectra with desired peaks
        if type(mz_list) is float or type(mz_list) is int:
            mz_list = [mz_list]
        if checks < 1:
            checks = 1
        else:
            pass
        if checks > len(mz_list):
            checks = len(mz_list)
        else:
            pass
        entries = 0
        neg_entries = 0
        new_mgf = []
        new_mgf_neg = []
        for entry in self.spectra:
            if entry.check_peaks(mz_list = mz_list, intensity_tol = intensity_tol, wobble_tol = wobble_tol) >= checks:
                new_mgf.append(entry)
                entries += 1
            else:
                new_mgf_neg.append(entry)
                neg_entries += 1
        self.pos_spectra = MGF()
        self.pos_spectra.spectra = new_mgf
        self.pos_spectra.filename = self.filename[:-4] + "_pos.mgf"
        self.neg_spectra = MGF()
        self.neg_spectra.spectra = new_mgf_neg
        self.pos_spectra.filename = self.filename[:-4] + "_neg.mgf"
        ###DEBUG TRACKING
        print(len(new_mgf))
        print(len(new_mgf)/len(self.spectra))
        self.filter_peaks = mz_list
        if negative == False:
            return new_mgf
        else:
            return new_mgf_neg

    def summary(self):
        string = ""
        string += "Spectra: {spectra}\n".format(spectra = len(self.spectra))
        av_pep = sum(float(x.pepmass) for x in self.spectra)/len(self.spectra) 
        string += "Average Precursor Mass: {av}\n\n".format(av = av_pep)
        if self.filter_peaks is None:
            string += "Not yet filtered\n\n"
        else:
            string += "Filter: {filter}\n".format(filter = ", ".join(map(str, self.filter_peaks)))
            #string += Tolerance/Checks:
            string += "Positive Spectra: {pos_spec}\n".format(pos_spec = len(self.pos_spectra.spectra))
            string += "Negative Spectra: {neg_spec}\n".format(neg_spec = len(self.neg_spectra.spectra))
            string += "Ratio: {pos_percent:.2f}% positive\n".format(pos_percent = len(self.pos_spectra.spectra)/(len(self.neg_spectra.spectra) + len(self.pos_spectra.spectra)) * 100)
            p_m = sum(float(x.pepmass) for x in self.pos_spectra.spectra)/len(self.pos_spectra.spectra) 
            string += "Average Precursor Mass (Positive Spectra): {pos_mass:.2f}\n".format(pos_mass = p_m)
            n_m = sum(float(x.pepmass) for x in self.neg_spectra.spectra)/len(self.neg_spectra.spectra) 
            string += "Average Precursor Mass (Negative Spectra): {neg_mass:.2f}\n\n".format(neg_mass = n_m)
        if len(self.unique_diffs) != 0: 
            string += "Mass Differences: {diffs}\n".format(diffs = len(self.unique_diffs))
        else:
            string += "Differences not calculated"
        return string

    def diff(self, diff_tol = 0.001):
        return list_diff(self, diff_tol = diff_tol)
    
    def list_diff(self, diff_tol = 0.001, intensity_tol = 0.2, filtered = False):
        """ calculates all mass differences within each spectra in the file

        Keyword Arguments:
        mgf_dict -- dictionary generated from a .mgf file
        diff_tol -- the difference in mz for two peaks to be considered different
        """
        #TODO: Frequency cutoff
        #TODO: Max spectra intensity threshold
        #iterate through mgf finding all differences and inserting them into diff_list
        #calculate mass differences
        print("start")
        print(datetime.now())
        for entry in self.spectra:
            valid_peaks = []
            spectra = entry.spectra
            int_thresh = max(spectra, key = lambda x:float(x[1]))[1]
            for pair in spectra:
                if pair[1] < int_thresh:
                    valid_peaks.append(pair[0])
                else:
                    pass                
            for row in range(len(valid_peaks) - 1):
                mass = valid_peaks[row]
                for column in range(row + 1, len(valid_peaks)):
                    self.diff_list.append([float(valid_peaks[column]) - float(mass), entry])
        diff_list = self.diff_list
        print(diff_list[0])
        print("diffs calced")
        print(datetime.now())
        #sort diff_list
        diff_list.sort(key = lambda x: x[0])
        print("sorted")
        print(datetime.now())
        print(len(diff_list))
        #iterate through diff list and make a count of all diffs
        diff_counts = []
        unique_diffs = []
        spectra_total = []
        spectra_list = []
        current = diff_list[0][0]
        count = 0
        max_mass = 0
        for mass in diff_list:
            if mass[0] > max_mass:
                max_mass = mass[0]
            else:
                pass
            if abs(mass[0] - current) <= diff_tol:
                count += 1
                spectra_list.append(mass[1])
            else:
                diff_counts.append(count)
                unique_diffs.append(current)
                spectra_total.append(spectra_list)
                spectra_list = [entry]
                current = mass[0]
                count = 1
        print("counted")
        self.diff_counts = diff_counts
        self.unique_diffs = unique_diffs
        self.diff_spectra = spectra_total
        #filter by fraction of max peak
        print(datetime.now())
        #return [diff_counts, unique_diffs]

    def tree_diff(self, diff_tol = 0.001):
        """ calculates all mass differences within each spectra in the file

        Keyword Arguments:
        mgf_dict -- dictionary generated from a .mgf file
        diff_tol -- the difference in mz for two peaks to be considered different
        """
        #TODO: Frequency cutoff
        #TODO: rounding/difference filter
        #TODO: add intensity filter
        #calculate mass differences
        start = datetime.now()
        for entry in self.spectra:
            spectra = entry.spectra
            for row in range(len(spectra) - 1):
                mass = spectra[row][0]
                for column in range(row + 1, len(spectra)):
                    difference = float(spectra[column][0]) - float(mass)
                    if self.diff_root is None:
                        self.diff_root = DifferenceNode(difference, 1, spectra, self)
                        self.diff_list.append(difference)
                    else:
                        self.diff_root.add_child(difference, spectra, diff_tol)
        print("start = ", start, "end = ", datetime.now())

    def find_diff(self, diff, tol = 0.001):
        found = []
        for current in range(len(self.unique_diffs)):
            if abs(self.unique_diffs[current] - diff) < tol:
                found.append([self.unique_diffs[current], self.diff_counts[current], self.diff_spectra[current]])
            else:
                pass
        return found

    def find_diff_tree(self, diff, tol = 0.001):
        if self.diff_root is not None:
            return self.diff_root.find(diff, tol)
        else:
            print("Differences not calculated. Calculate Differences using the \"diff\" method.")

class Spectra:
    def __init__(self, filename, title, pepmass, charge, scans, rt, spectra, rtsecs = True):
        self.filename = filename
        self.title = title
        self.pepmass = pepmass
        self.charge = charge
        self.scans = scans
        self.rtinseconds = rt   
        self.spectra = spectra
        #track whether rt is in seconds or minutes
        self.rtsecs = rtsecs 

    def __str__(self):
        #return whatever i want it to print
        #construct string
        return "info stuff"

    def show(self):
        print("File: " + self.filename)
        print("Title: " + self.title)
        print("Pepmass: " + self.pepmass)
        print("Charge: " + self.charge)
        print("Scans: " + self.scans)
        print("RT in Seconds: " + self.rtinseconds)
        print("Spectra: ")
        for i in self.spectra:
            print(i[0] + " , " + i[1])

    def check_peaks(self, mz_list, intensity_tol, wobble_tol):
        """Checks whether a given spectra contains peaks of specified m/z values
        
        Keyword Arguments:
        mz_list -- a list of desired m/z values to check spectra for
        intensity_tol -- intensity relative to the max peak in each spectra that a peak must be above to be detected
        wobble_tol -- value above and below the specified m/z values peaks can be detected
        """
        if type(mz_list) is float or type(mz_list) is int:
            mz_list = [mz_list]
        #Check if list?
        mz_presence = []
        mz_occur = 0
        max_peak = 0
        for mz in mz_list:
            mz_peak = 0
            for peak in self.spectra:
                #check for max
                if float(peak[1]) > max_peak:
                    max_peak = float(peak[1])
                else:
                    pass
                #Tolerance for peak wobble?
                if float(peak[0]) - mz < wobble_tol and float(peak[0]) - mz >= 0 and float(peak[1]) > mz_peak:
                    mz_peak = float(peak[1])
                else:
                    pass
            if mz_peak > max_peak*intensity_tol:
                mz_presence.append(True)
                mz_occur += 1
            else:
                mz_presence.append(False)
            return mz_occur

        def max_intensity(self):
            return max(self.spectra[1].spectra, key = lambda x:float(x[1]))[1]

class DifferenceNode:
    def __init__(self, difference, count, spectra, mgf = None):
        self.difference = difference
        self.count = count
        self.spectra = spectra #???
        self.mgf = mgf

        self.left_child = None
        self.right_child = None

    def __str__(self):
        string = ""
        string += "Difference: {difference:.3f}\n".format(difference = self.difference)
        string += "Count: {count}\n".format(count = self.count)
        return string

    def add_child(self, new_diff, spectra, diff_tol = 0.001):
        #TODO: add tolerance
        if diff_tol >= abs(self.difference - new_diff):
            self.count += 1
        elif new_diff > self.difference:
            if self.left_child == None:
                self.left_child = DifferenceNode(new_diff, 1, self.spectra, self.mgf)
            else:
                self.left_child.add_child(new_diff, spectra, diff_tol)
                self.mgf.diff_count += 1
                #self.mgf.diff_list.append(self.left_child)
        elif new_diff < self.difference:
            if self.right_child == None:
                self.right_child = DifferenceNode(new_diff, 1, self.spectra, self.mgf)
            else:
                self.right_child.add_child(new_diff, spectra, diff_tol)
                self.mgf.diff_count += 1
                #self.mgf.diff_list.append(self.right_child)

    def find(self, difference, tolerance = 0.001):
        #TODO: add tolerance
        if tolerance >= abs(difference - self.difference):
            return self
        elif difference > self.difference:
            if self.left_child is None:
                print("Not Found")
            else:
                return self.left_child.find(difference)
        elif difference < self.difference:
            if self.right_child is None:
                print("Not Found")
            else:
                return self.right_child.find(difference)
###############################################################
def extract_mgf(mgf):
    #TODO: handle - mass charges
    """Read a .mgf into a dictionary
    
    Keyword Arguments:
    mgf -- the filepath of the .mgf file to be extracted
    """
    mgf_file = open(mgf,  "r")
    mgf_split = str(mgf_file.read()).split("END IONS")
    mgf_list = []
    for i in range(len(mgf_split) - 1):
        mgf_part = mgf_split[i].split("\n")
        if i != 0:
            mgf_part = mgf_part[2:]
        spectra = mgf_part[6:-1]
        for k in range(len(spectra)):
            spectra[k] = spectra[k].split(" ")
        mgf_dict = {
            "FILE": mgf,
            "TITLE": mgf_part[1][12:],
            "PEPMASS": mgf_part[2][8:],
            "CHARGE": mgf_part[3][7:-1],
            "SCANS": mgf_part[4][6:],
            "RTINSECONDS": mgf_part[5][12:],
            "SPECTRA": spectra
        }
        mgf_list.append(mgf_dict)
    return mgf_list

def extract_peaks(mgf_list, mz_list, intensity_tol = 0.2, wobble_tol = 0.01, checks = 1):
    """Takes a .mgf file converted into a dictionary by extract_mgf() and extracts all spectra with peaks of      a desired m/z value and 
    intensity relative to the max intensity in the spectra
    
    Keyword Arguments:
    mgf_list -- dictionary generated from a .mgf file
    mz_list -- a list of desired m/z values to check spectra for
    intensity_tol -- intensity relative to the max peak in each spectra that a peak must be above to be           detected
    wobble_tol -- value above and below the specified m/z values peaks can be detected
    checks -- the number of peaks that a spectra must have with the specified m/z value to be detected, a value below 1 will require 1 or more peaks to be present 
                and a value greater than the amount of specified peaks will require all peaks to be present
    """
    #Extract spectra with desired peaks
    if checks < 1:
        checks = 1
    else:
        pass
    if checks > len(mgf_list):
        checks = len(mgf_list)
    else:
        pass
    entries = 0
    neg_entries = 0
    new_mgf = []
    new_mgf_neg = []
    for entry in mgf_list:
        if find_peaks(mz_list = mz_list, spectra = entry["SPECTRA"], intensity_tol = intensity_tol, wobble_tol = wobble_tol) >= checks:
            new_mgf.append(entry)
            entries += 1
        else:
            new_mgf_neg.append(entry)
            neg_entries += 1
    ###DEBUG TRACKING
    print(len(new_mgf))
    print(len(new_mgf)/len(mgf_list))
    with open("Subset Stats.txt", "a") as f:
        f.write(str(mz_list) + "\n")
        f.write(str(len(new_mgf)) + "/" + str(len(mgf_list)) + "\n")
        f.write(str(len(new_mgf)/len(mgf_list)) + "\n")
        f.write("\n")
    return new_mgf, new_mgf_neg #, stats?

def find_peaks(mz_list, spectra, intensity_tol = 0.2, wobble_tol = 0.01):
    """Checks whether a given spectra contains peaks of specified m/z values
    
    Keyword Arguments:
    mz_list -- a list of desired m/z values to check spectra for
    spectra -- spectra to be searched presented as a list of lists
    intensity_tol -- intensity relative to the max peak in each spectra that a peak must be above to be detected
    wobble_tol -- value above and below the specified m/z values peaks can be detected
    """
    if type(mz_list) is float or type(mz_list) is int:
        mz_list = [mz_list]
    #Check if list?
    mz_presence = []
    mz_occur = 0
    max_peak = 0
    for mz in mz_list:
        mz_peak = 0
        for peak in spectra:
            #check for max
            if float(peak[1]) > max_peak:
                max_peak = float(peak[1])
            else:
                pass
            #Tolerance for peak wobble?
            if float(peak[0]) - mz < wobble_tol and float(peak[0]) - mz >= 0 and float(peak[1]) > mz_peak:
                mz_peak = float(peak[1])
            else:
                pass
        if mz_peak > max_peak*intensity_tol:
            mz_presence.append(True)
            mz_occur += 1
        else:
            mz_presence.append(False)
        return mz_occur
    
def save_mgf(mgf, filename):
    """saves a dictionary gernerated from a.mgf file back into a .mgf file

    Keyword Arguments:
    mgf -- dictionary generated from a .mgf file
    filename -- filename to give the new .mgf file
    """
    filename_mgf = str(filename) + ".mgf"
    #Write to new mgf file
    #Blank file
    with open(filename_mgf, "w") as f:
        f.write("")
    new_file_write = open("filename_mgf", "w")
    new_file_write.write("")
    #Append each line to file
    for entry in mgf:
        with open(filename_mgf, "ab") as f:
            f.write(bytearray(b"BEGIN IONS\n"))
            f.write(bytearray("TITLE=index=" + entry["TITLE"] + "\n", "utf-8"))
            f.write(bytearray("PEPMASS=" + entry["PEPMASS"] + "\n", "utf-8"))
            f.write(bytearray("CHARGE=" + entry["CHARGE"] + "+\n", "utf-8"))
            f.write(bytearray("SCANS=" + entry["SCANS"] + "\n", "utf-8"))
            f.write(bytearray("RTINSECONDS=" + entry["RTINSECONDS"] + "\n", "utf-8"))
            for peak in entry["SPECTRA"]:
                f.write(bytearray(str(peak[0]) + " " + str(peak[1]) + "\n", "utf-8"))
            f.write(b"END IONS\n\n")
    
def extract_to_mgf(mgf_file, mz_list, intensity_tol = 0.2, wobble_tol = 0.01, new_filename = "filter_result", checks = 1):
    """
    """
    mgf_list = extract_mgf(mgf_file)
    new_mgf, new_mgf_neg = extract_peaks(mgf_list, mz_list, intensity_tol, wobble_tol, checks)
    save_mgf(new_mgf, new_filename)
    save_mgf(new_mgf_neg, new_filename + "_neg")

def all_diffs(mgf_dict, diff_tol = 0.001):
    """ calculates all mass differences within each spectra in the file

    Keyword Arguments:
    mgf_dict -- dictionary generated from a .mgf file
    diff_tol -- the difference in mz for two peaks to be considered different
    """
    #TODO: Frequency cutoff
    #RobParkerSpike.mgf - ~ 4 minutes
    diff_list = []
    #iterate through mgf finding all differences and inserting them into diff_list
    diff_list = []
    #calculate mass differences
    for entry in mgf_dict:
        spectra = entry["SPECTRA"]
        for row in range(len(spectra) - 1):
            mass = spectra[row][0]
            for column in range(row + 1, len(spectra)):
                diff_list.append(float(spectra[column][0]) - float(mass))
    print("diffs calced")
    print(datetime.now())
    #sort diff_list
    diff_list.sort()
    print("sorted")
    print(datetime.now())
    print(len(diff_list))
    #iterate through diff list and make a count of all diffs
    diff_counts = []
    unique_diffs = []
    current = diff_list[0]
    count = 0
    max_mass = 0
    for mass in diff_list:
        if mass > max_mass:
            max_mass = mass
        else:
            pass
        if abs(mass - current) <= diff_tol:
            count += 1
        else:
            diff_counts.append(count)
            unique_diffs.append(current)
            current = mass
            count = 1
    print("counted")
    #filter by fraction of max peak
    print(datetime.now())
    return diff_counts, unique_diffs

def all_diffs_collate(mgf_dict, diff_tol = 0.0001):
    """ calculates all mass differences within each spectra in the file

    Keyword Arguments:
    mgf_dict -- dictionary generated from a .mgf file
    diff_tol -- the difference in mz for two peaks to be considered different
    """
    #TODO: Frequency cutoff
    #RobParkerSpike.mgf - ~ 4 minutes
    diff_list = []
    #iterate through mgf finding all differences and inserting them into diff_list
    diff_list = []
    #calculate mass differences
    for entry in mgf_dict:
        spectra = entry["SPECTRA"]
        for row in range(len(spectra) - 1):
            mass = spectra[row][0]
            for column in range(row + 1, len(spectra)):
                diff_list.append(float(spectra[column][0]) - float(mass))
    print("diffs calced")
    print(datetime.now())
    #sort diff_list
    diff_list.sort()
    print("sorted")
    print(datetime.now())
    print(len(diff_list))
    #iterate through diff list and make a count of all diffs
    diff_counts = []
    unique_diffs = []
    current = diff_list[0]
    count = 0
    max_mass = 0
    for mass in diff_list:
        if mass > max_mass:
            max_mass = mass
        else:
            pass
        if abs(mass - current) <= diff_tol:
            count += 1
        else:
            diff_counts.append(count)
            unique_diffs.append(current)
            current = mass
            count = 1
    ######################################
    #COLLATE
    #set initial value (difference_list = 0)
    #for difference in difference_list
        #check if local maximum
            #if yes make local maximum
                #Check distance of current from local maximum
                #If distance > tolerance, collate range(initial, current - 1) into a single peak
                #else continue
            #else 
                #check distance of initial from local maximum
                    #???
                    #If distance > tolerance, collate range(initial, current - 1) into a single peak
                    #else continue
    ######################################
    print("counted")
    #filter by fraction of max peak
    print(datetime.now())
    return diff_counts, unique_diffs

def filter_diffs(the_diffs, did_diffs, tolerance = 0.01):
    #filter differences
    print(len(did_diffs))
    filter_the = []
    filter_did = []
    threshold = max(did_diffs)*tolerance
    pos_count = 0
    neg_count = 0
    diff = 0
    while diff < len(the_diffs):
        if did_diffs[diff] < threshold:
            neg_count += 1
        else:
            filter_the.append(the_diffs[diff])
            filter_did.append(did_diffs[diff])
            pos_count += 1
    print(neg_count)
    return filter_the