\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{assignment}[2016/02/22]

\LoadClass[10pt, twoside]{article}

\usepackage{fancyhdr,rotating,amsmath}
\usepackage{graphicx}
\usepackage{listings}
\usepackage{anysize}
\usepackage{array}

\marginsize{2.0cm}{2.0cm}{2.0cm}{3.0cm}
%\pagenumbering{gobble}

\setlength{\parindent}{0ex}
\setlength{\parskip}{1em}

% Title
\newcommand{\mytitle}[1]{\gdef\@mytitle{#1}}

\pagestyle{fancy}
 \let\headrule=\empty
 \let\footrule=\empty
\fancyhead{}
\fancyhead[CO]{ {\Large \bfseries \@mytitle}}

\lstset{language=Python}

\newcommand{\code}{\texttt}

\endinput
