import pathlib
import tkinter as tk
from tkinter import DISABLED, END, filedialog
from tkinter.font import NORMAL
from tkinter.ttk import Combobox
import ube3_apa
import threading
import sys

class PrintConsole():
    def __init__(self, box):
        self.box = box

    def write(self, text):
        self.box.insert(tk.END, text)

    def flush(self):
        pass

if __name__ == '__main__':

    window = tk.Tk()
    window.title("UbE3-APA")
    window.geometry('600x600')

    lblInput = tk.Label(window, text="Select Input File for Site Ratio (required)")
    lblInput.grid(column=0, row=0)
    txtInput = tk.Entry(window,width=50)
    txtInput.grid(column=1, row=0)
    
    
    
    def pronormyes():
            btnBrowseInput_pro['state'] = tk.NORMAL
    def pronormno():
            btnBrowseInput_pro['state'] = tk.DISABLED
            
            
    lblpronorm = tk.Label(window, text="Normalize by protein ratio")
    lblpronorm.grid(column=0, row=1)
    pronorm = tk.IntVar()
    pronorm.set(0)
    rpronormTrue = tk.Radiobutton(window, text="Yes", variable=pronorm, value=1, command = pronormyes)
    rpronormFalse = tk.Radiobutton(window, text="No", variable=pronorm, value=0, command = pronormno)
    rpronormTrue.grid(column=1, row=1)
    rpronormFalse.grid(column=2, row=1)
    
    
    # lblpronorm = tk.Label(window, text="Normalize by protein ratio")
    # lblpronorm.grid(column=0, row=1)
    # buttonnormyes = tk.Button(window, text="Yes", command = switchnormState)
    # buttonnormyes.grid(column=1, row=1)
    # buttonnormno = tk.Button(window, text="No", command = switchnormState, state=tk.DISABLED)
    # buttonnormno.grid(column=2, row=1)


    
    lblInput = tk.Label(window, text="Select Input File for Protein Ratio (optional)")
    lblInput.grid(column=0, row=2)
    txtInput_pro = tk.Entry(window,width=50)
    txtInput_pro.grid(column=1, row=2)

    lblOutput = tk.Label(window, text="Set Output Directory (required)")
    lblOutput.grid(column=0, row=3)
    txtOutput = tk.Entry(window,width=50)
    txtOutput.grid(column=1, row=3)
    
    def getInput():
        #inputDirectory = filedialog.askdirectory(initialdir='/', title="Select Input Directory")
        inputFile = filedialog.askopenfile(mode='r', filetypes=[('TXT Files', '*.txt')])
        inputFilePath = pathlib.Path(inputFile.name)
        txtInput.delete(0,END)
        txtInput.insert(0, inputFilePath)
    def getInput_pro():
        #inputDirectory = filedialog.askdirectory(initialdir='/', title="Select Input Directory")
        inputFile = filedialog.askopenfile(mode='r', filetypes=[('TXT Files', '*.txt')])
        inputFilePath = pathlib.Path(inputFile.name)
        txtInput_pro.delete(0,END)
        txtInput_pro.insert(0, inputFilePath)
    def setOutput():
        outputDirectory = filedialog.askdirectory(initialdir='/', title="Select Output Directory")
        txtOutput.delete(0,END)
        txtOutput.insert(0, outputDirectory)

    btnBrowseInput = tk.Button(window, text="Browse", command=getInput)
    btnBrowseInput.grid(column=2, row=0)
    
    btnBrowseInput_pro = tk.Button(window, text="Browse", command=getInput_pro,state=tk.DISABLED)
    btnBrowseInput_pro.grid(column=2, row=2)

    btnBrowseOutput = tk.Button(window, text="Browse", command=setOutput)
    btnBrowseOutput.grid(column=2, row=3)


  
    
    lblOutput = tk.Label(window, text="Experiment Label in MaxQuant Results")
    lblOutput.grid(column=0, row=4)
    txtExplabel = tk.Entry(window,width=50)
    txtExplabel.grid(column=1, row=4)



    # lblOutput = tk.Label(window, text="Input Type")
    # lblOutput.grid(column=0, row=5)
    # inputTypeComboBox = Combobox(window, values=("UniprotAC", "protein", "gene symbol", "gene"))
    # inputTypeComboBox.grid(column=1, row=5)
    inputType = 'UniprotAC'


    lblGrouped = tk.Label(window, text="Grouped Mode")
    lblGrouped.grid(column=0, row=6)
    grouped = tk.IntVar()
    grouped.set(1)
    rGroupTrue = tk.Radiobutton(window, text="True", variable=grouped, value=1)
    rGroupFalse = tk.Radiobutton(window, text="False", variable=grouped, value=0)
    rGroupTrue.grid(column=1, row=6)
    rGroupFalse.grid(column=2, row=6)



    lblRatio = tk.Label(window, text="Export Ratio Tables")
    lblRatio.grid(column=0, row=7)
    outputRatio = tk.IntVar()
    outputRatio.set(1)
    rOutputRatioTrue = tk.Radiobutton(window, text="True", variable=outputRatio, value=1)
    rOutputRatioFalse = tk.Radiobutton(window, text="False", variable=outputRatio, value=0)
    rOutputRatioTrue.grid(column=1, row=7)
    rOutputRatioFalse.grid(column=2, row=7)



    lblLog2Trans = tk.Label(window, text="Log2Trans Iuput data")
    lblLog2Trans.grid(column=0, row=8)
    outputLog2Trans = tk.IntVar()
    outputLog2Trans.set(1)
    rLog2TransTrue = tk.Radiobutton(window, text="True", variable=outputLog2Trans, value=1)
    rLog2TransFalse = tk.Radiobutton(window, text="False", variable=outputLog2Trans, value=0)
    rLog2TransTrue.grid(column=1, row=8)
    rLog2TransFalse.grid(column=2, row=8)

    consoleOutput = tk.Text(window, height=20, width=50)

    def switchState(b):
        if b['state'] == NORMAL:
            b['state'] = DISABLED
        else:
            b['state'] = NORMAL

    def run(inputDirectory, inputType, outputDirectory, exp_label, pronorm, inputDirectory_pro, grouped, ratio, log2trans, t, b):
        console = PrintConsole(t)
        sys.stdout = console
        def boolify(n):
            n = n.get()
            if n == 0:
                return False
            else:
                return True
        def threadWorker(inputDirectory, inputType, outputDirectory, exp_label, pronorm, inputDirectory_pro, grouped, ratio, log2trans, b):
            print('Analysis Started. Please be patient...')
            if(boolify(pronorm) == True):
                paths = ube3_apa.maxqt2ube3( output_dir=pathlib.Path(outputDirectory), site_ratio_dir=pathlib.Path(inputDirectory), pro_ratio_dir=pathlib.Path(inputDirectory_pro), exp_label=exp_label)
                
            else:
                paths = ube3_apa.maxqt2ube3( output_dir=pathlib.Path(outputDirectory), site_ratio_dir=pathlib.Path(inputDirectory), pro_ratio_dir='', exp_label=exp_label)
            exp_label = exp_label.replace("/", "")
            ube3_apa.e3enrich(paths[0],inputType, pathlib.Path(outputDirectory), exp_label=exp_label, proratio_dir=paths[1], grouped=boolify(grouped), ratio_output=boolify(ratio), log2trans=boolify(log2trans)) 
            print('Analysis Completed')
            switchState(b)
        thread = threading.Thread(target = lambda:threadWorker(inputDirectory, inputType, outputDirectory, exp_label, pronorm, inputDirectory_pro, grouped, ratio, log2trans, b))
        thread.start()
    
    btnSubmit = tk.Button(window, text="Submit", command=lambda:[switchState(btnSubmit), run(txtInput.get(), inputType, txtOutput.get(), txtExplabel.get(), pronorm, txtInput_pro.get(), grouped, outputRatio, outputLog2Trans, consoleOutput, btnSubmit)])
    btnSubmit.grid(column=2, row=9)
    lblConsole = tk.Label(window, text="Console")
    lblConsole.grid(column=0, row=10)
    consoleOutput.grid(column = 0, columnspan=4, row=11)

    window.mainloop()