shshs
# 🚀 Getting Started

## 1. Check Python Installation

Open a terminal and run:

```bash
python --version
```

or

```bash
python3 --version
```

You should see Python 3.10 or newer.

If Python is not installed:

- Windows: https://www.python.org/downloads/
- macOS: https://www.python.org/downloads/
- Linux: Use your package manager or install from python.org

---

## 2. Clone the Repository

```bash
git clone git@gitlab.stfc.ac.uk:isis-accelerator-controls/isis_wx_26.git
cd gui/tutorials/streamlit_tutorial/
```

Alternatively, download and extract the ZIP file if you do not have Git installed.

---

## 3. Create a Virtual Environment (Recommended)

A virtual environment prevents package conflicts between projects.

### Windows

```bash
python -m venv venv
venv\Scripts\activate
```

### Linux / macOS

```bash
python3 -m venv venv
source venv/bin/activate
```

You should now see `(venv)` at the beginning of your terminal prompt.

---

## 4. Upgrade pip

Before installing packages:

```bash
python -m pip install --upgrade pip
```

---

## 5. Install Project Requirements

Install all required dependencies:

```bash
pip install -r requirements.txt
```

This installs:

- Streamlit
- Pandas
- NumPy
- Plotly
- Seaborn
- Matplotlib

---

## 6. Verify the Installation

Run:

```bash
python -c "import streamlit, pandas, numpy, plotly, seaborn, matplotlib; print('Installation successful!')"
```

If no errors appear, you're ready to go.

---

## 7. Run the Application

Start the Streamlit server:

```bash
streamlit run tutorial.py
```

After a few seconds, Streamlit will automatically open a browser window.

If it doesn't open automatically, visit:

```text
http://localhost:8501
```

---

# Common Issues

## Command Not Found: streamlit

Try:

```bash
python -m streamlit run app.py
```

---

## ModuleNotFoundError

The required packages are not installed.

Run:

```bash
pip install -r requirements.txt
```

---

## Permission Errors

Try:

```bash
python -m pip install --upgrade pip
```

and make sure your virtual environment is activated.

---

## Virtual Environment Not Activating

### Windows PowerShell

If you receive an execution policy error:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

Then try:

```powershell
venv\Scripts\activate
```

again.

---

# Project Structure

```text
streamlit_tutorial/
│
├── app.py
├── requirements.txt
└── README.md
```

---

# Closing the Virtual Environment

When you are finished working:

```bash
deactivate
```

To return later:

### Windows

```bash
venv\Scripts\activate
```

### Linux / macOS

```bash
source venv/bin/activate
```

# Foot Note for the EPICS ARCHIVER TEAM

In the browser you can put the following in the search bar:

## Request data for a specific PV
```bash
http://athena.isis.rl.ac.uk:9506/data?pv=HELREC::LIQUIFIER:PT112:READ&from=2025-06-19T16:40:00.000000Z&to=2025-06-19T17:10:00.000000Z
```
Break down of the request:
* http://athena.isis.rl.ac.uk:9506/ - IP for the API endpoints
* pv=HELREC::LIQUIFIER:PT112:READ - PV stands for Process Variable **HELREC::LIQUIFIER:PT112:READ** is the PV we are getting data from
* from=2025-06-19T16:40:00.000000Z - all the data will be older then this timestamp, some times there is no data at the initial timestamp as such this does not always represent the timestamp of the first element.
* to=2025-06-19T17:10:00. - it aims to get data up to this timestamp

## PV Status
Use the following to check if a PV exists
```bash
http://athena.isis.rl.ac.uk:9506/getPVStatus?pv=HELREC::LIQUIFIER:PT112:READ
```
```bash
http://athena.isis.rl.ac.uk:9506/getPVStatus?pv=IDONTEXIST:LOL
```

## Regular Expresion Search
Use the following to find all the relevant PVs we have in the system
```bash
http://athena.isis.rl.ac.uk:9506/glob?pv=HELREC::LIQUIFIER:FT*
```
