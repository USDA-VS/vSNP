---
layout: "default"
---

<h1><p style="text-align: center">Windows 10</p></h1>

-----
<br>

Windows 10 Ubuntu app setup
=================

Install Ubuntu 18 LTS from Microsoft app store

Launch app

If needed `Enable-WindowsOptionalFeature`.  Follow instructions provided by link.  Run Powershell as administrator (right-click on Powershell).

Launch Ubuntu app.  You will be prompted to make a username and password.

Your home directory will be: ~$ /home/\<username\>

Download Linux Anaconda Installer, `https://www.anaconda.com/download/#linux`, Right click on "64-Bit(86) Installer" link and "Copy link address".

`~$ wget <copied link>`

`~$ bash ./Anaconda3*` and follow prompts (installation will take a while, if installation stalls, press \<enter\>)

Conda should prepend this or something similar to your ~/.bashrc:

```. ${HOME}/anaconda3/etc/profile.d/conda.sh`
conda activate```

Follow [Linux/Mac Setup](https://usda-vs.github.io/vSNP/setup.html). Start by cloning vSNP repository.

A note about the subsystem...

From within the Windows 10 Ubuntu subsystem terminal the directory structure is the familiar Linux structure.  Directories under root are as expected as with any Linux environment:  `/home`, `/usr`, `etc`.  The Windows directory structure is under `/mnt`.  The Windows c drive is a subdirectory under `/mnt`.  Therefore your Windows 10 Desktop may have a path such as: `/mnt/c/Users/<windows username>/OneDrive/Desktop`
