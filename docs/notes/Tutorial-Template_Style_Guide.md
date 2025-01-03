<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# The NRPy+ Jupyter Tutorial Style Guide / Template
## Authors: Brandon Clark, Zach Etienne, & First Last
### Formatting improvements courtesy Brandon Clark

<font color='red'>**This is a warning message, in red text and bolded, to warn anyone using the module that it is, for example, actively in development, not yet validated, etc. Warning messages are optional.**</font>

## This module implements a template designed by Brandon Clark to be used as a style guide for all tutorial notebooks within NRPy+.

### Items in Markdown code contained within "</>" are not included within the output (double click this box to see what I mean). To the run Markdown code simply hit "Shift + Enter" or the "Run" button above. 

<font color='green'>**This text discusses how a module has been validated against other existing code or modules. This text is given a green font color and is bolded. See how to bold and make text different colors in the Markdown code.**</font>

### </list_source_code> NRPy+ Source Code for this module:
1. [Template_Style_Guide.py](../edit/Template_Style_Guide.py); [\[**tutorial**\]](Tutorial-Template_Style_Guide.ipynb) </description_here> This is where you would describe what purpose this source code serves in this module. Read how to correctly link to these source code files/tutorial notebooks later in []. 
1. </additional_source_code_links__here>

## Introduction:
Here you write an introduction that discusses in slight detail the framework of this tutorial notebook. Here you may reference external works or websites on which pieces of your module rely. It is often helpful to include an enumerated algorithm to highlight this module's processes. Within the algorithm, you may refer to where source code is implemented as a part of this module. </optional>

The entire </made_up>algorithm is outlined below, with NRPy+-based components highlighted in <font color='green'>green</font>.

1. Constructing a Table of Contents
1. 1. Discussing [Markdown Linking Protocol](https://medium.com/@sambozek/ipython-er-jupyter-table-of-contents-69bb72cf39d3)
    1. Linking to sections internally within the module 
    1. Linking to external sources
1. No parts of this template tutorial notebook rely on <font color='green'>NRPy+-based components</font>
1. Converting Jupyter Notebook to output LaTex PDF

</optional>
You could also write your introduction to include subsections preceded by ###. 

### Introduction subsection:
Include information relevant to this subsection here.

## </other>  Other (Optional): 
You may include any number of items here within the first box of the tutorial notebook, but I suggest being minimalistic when you can. Other sections that have been included in other tutorial modules are as follows

### Note on Notation:
When using a new type of notation for the first time within the NRPy+ tutorial, you may want to include some notes on that here.

### Citations:
This is a great place to list out the references you link to within the module with actual citations. 

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows.

0. [Preliminaries](#prelim): This is an optional section
1. [Step 1](#linking): The Markdown Linking Protocol </header_section>
    1. [Step 1.a](#internal_links) Internal linking with the Jupyter notebook, Table of Contents </subsection>
    1. [Step 1.b](#external_links): External linking outside of Jupyter notebook </subsection>
        1. [Step 1.b.i](#nrpy_links): Linking to other files/modules within NRPy+ </subsubsection>
1. [Step 2](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file


The Table of Contents (ToC) plays a significant role in the formatting of your module. The above ToC is for this module, but I have constructed it in a way such that you should see all of the important details for any module you need to write. If you choose to include a preliminaries section, enumerate it with the "0." All other sections, subsections, and sub-subsections can be enumerated with the 1. Jupyter/LaTex will handle their own numbering/lettering scheme. It is important when creating subsections and sub-subsections that you indent as seen in the Markdown code. The text colors vary for the level section you're assigning within the Markdown code. When writing within the brackets to specify a step number, the following scheme is to be used.

* Header Sections: Step 1, Step 2, Step 3
* Subsections: Step 1.a, Step 1.b, Step 1.c
* Sub-subsections: Step 1.a.i, 1.a.ii, 1.a.iii </roman_numerals>

If for some reason you go more then three levels deep in your sectioning, I would suggest finding a way to reorganize your sectioning to prevent that, or ask Zach Etienne what the next level of labeling for Steps should be. We will talk about the other components within the Markdown Code for the ToC in [Step 1.a](#internal_links). The only text within the ToC section of this module should be the ToC code itself and what precedes it.

I also suggest that the titles for the Steps you include here following the ":" match the titles you use throughout your module.

<a id='prelim'></a>

# Preliminaries: This is an optional section \[Back to [top](#toc)\]
$$\label{prelim}$$ 

This section is a great chance to include textual verbiage that might have been too specific for the introduction section but serves as a beneficial setup for the remainder of the module. For instance, you may want to define quantities here, express important equations, and so on. I suggest that the Preliminaries section not be followed by any Python code blocks and remains simply a block of information for users to refer back to.


<a id='linking'></a>

# Step 1: The Markdown Linking Protocol \[Back to [top](#toc)\]
$$\label{linking}$$

We have already within this template had to link to sources both internally within this module, externally to other components of the NRPy+ tutorial, as well as externally to additional web sources. The next few sections discuss how this is done. It is important to know that any linking is done by combining brackets and parentheses "\[ \]()" with the desired input in each. 

On another note, main sections like this have their titles preceded by a single #. As you will see, for every deeper layer of sectioning, an additional # is appended, reducing the size of the text.

<a id='internal_links'></a>

## Step 1.a: Internal linking with the Jupyter notebook, Table of Contents \[Back to [top](#toc)\]
$$\label{internal_links}$$

A great resource for how to construct a Table of Contents is:
https://medium.com/@sambozek/ipython-er-jupyter-table-of-contents-69bb72cf39d3.

The Table of Contents is a family of internal links. To link internally we first have to specify an ***anchor tag*** which is the text within the parentheses preceded by a # (See ToC Markdown code). For instance, the anchor tag for this subsection is `internal_links`. So, for a particular Step within the Table of Contents you specify the Step title in brackets (e.g., <font color='blue'>[Step 1.a]</font>), appended by the anchor tag in parentheses preceded by a # (e.g., <font color='red'>(#internal_links)</font>), followed by a ":" and the Step description (e.g., Internal linking with the Jupyter notebook, Table of Contents). Look at the Markdown code for the Table of Contents for a few examples. 

**Important Note**: The anchor tags cannot be anything that you want. Anchor tags must be entirely lowercase and contain no spaces. Numbers are fine as well as underscores, but not capitalization. I suggest making the anchor tags have significant meaning to the section there tied to, instead of making one that reads "step1a". The reason I say this is because if you ever need to resection your module, the tags won't all need to be changed as well if you give each one a unique name.  

All we have done so far is establish anchor tags and clickable links within the Table of Contents, but how do we establish the link to the specific section within the module? Opening up the Markdown code for this section you will see a line of code above the title and a line of code directly below the title. These are the answers to the question. Each section requires these components to be included for both the Jupyter Notebook and LaTex internal linking. Make sure the top line of the Markdown code has a space between it and the title. Similarly, the code directly beneath the title needs  space below it as well, separated from the main body of text (see above in Markdown code).

**Important Note**: Links do not work unless the two sections which are linked have been run.

The Table of Contents is now linked to this section, and you may have already noticed that this section, and all others, are linked back to the Table of Contents using the Markdown code in line at the end of the section title. This is exceedingly convenient for modules of great length. It may also be convenient when you're in a particular subsection and you wish to just return to the header section. This is accomplished using a bracket parentheses \[\]() pairing like so (see this in Markdown code). Go back to [Step 1](#linking).

Lastly, you would more often than not write a code block below implementing what was discussed in this section. This isn't always necessary, some header sections plainly serve as a setup for subsections that will contain all of the necessary coding components.  


```python
# This is the code block corresponding to Step 1.a: Internal linking within the Jupyter notebook, Table of Contents
print("We have successfully learned how to code internal links using Markdown Linking Protocol!!!")
```

    We have successfully learned how to code internal links using Markdown Linking Protocol!!!


<a id='external_links'></a>

## Step 1.b: External linking outside of this module \[Back to [top](#toc)\]
$$\label{external_links}$$

To link outside of this particular module we still use bracket parentheses \[ \]() pairings. Since the links are not internal, we no longer need the # symbol and anchor tags. Instead, you need an actual link. For instance, look at your Markdown code to see how we link this [website](https://medium.com/@sambozek/ipython-er-jupyter-table-of-contents-69bb72cf39d3) to a line of text. Of course, web links will simply work on their own as a hyperlink, but often you may need to link to multiple external sources and do not want all of the individual addresses clogging up the body of your text. 


```python
# This is the code block for Step 1.b: External linking outside of Jupyter notebook
print("Be efficient in how you link external sources, utilize []() pairs!!!")
```

    Be efficient in how you link external sources, utilize []() pairs!!!


<a id='nrpy_links'></a>

### Step 1.b.i: Linking to other files/modules within NRPy+ \[Back to [top](#toc)\]
$$\label{nrpy_links}$$

Other useful external sources we would like to link to are the existing files/modules within NRPy+. To do this we again resort to the \[ \]() pair. By simply typing the file name into the parentheses, you can connect to another [tutorial module](Tutorial-Template_Style_Guide.ipynb) (see Markdown). To access a .py file, you want to type the command ../edit/
followed by the file location. For instance, here is the [.py file](../edit/Template_Style_Guide.py) for this notebook (see Markdown). 



```python
# This is the code block for Step 1.b.i: Linking to other files/modules within NRPy+
print("Template_Style_Guide.py is an empty file...")
```

    Template_Style_Guide.py is an empty file...


<a id='latex_pdf_output'></a>

# Step 2: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-Template_Style_Guide.pdf](Tutorial-Template_Style_Guide.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)

**Important Note**: Make sure that the file name is right in all six locations, two here in the Markdown, and four in the code below. 

* Tutorial-Template_Style_Guide.pdf
* Tutorial-Template_Style_Guide.ipynb
* Tutorial-Template_Style_Guide.tex


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Template_Style_Guide")
```

    Created Tutorial-Template_Style_Guide.tex, and compiled LaTeX file to PDF
        file Tutorial-Template_Style_Guide.pdf

