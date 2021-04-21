Please follow the tutorial with images to transform your key, connect via Putty and use a scp client like winSCP

See here: https://github.com/clairemerot/physalia_adaptation_course/blob/2021/connect_AWS_windows.pdf


If you are using a Windows machine, you will need to log on using PuTTY since there is no native ssh client. PuTTY does not natively support the private key format (.pem) needed to login to our Amazon cloud instance. You first need to convert the private key that we gave to you to a key that PuTTY can read. PuTTY has a tool named PuTTYgen, which can convert keys to the required PuTTY format (.ppk). When you installed PuTTY, it will also have installed PuTTYgen.

First, start PuTTYgen (for example, from the Start menu, choose All Programs > PuTTY > PuTTYgen). Then select RSA and click on Load:

In the new window that pops up, Change “PuTTY Private Key Files” to “All Files” to allow you to find your pem file.

Then you will save your private key. Click on YES to dismiss the Warning.

Save your key with the same name as the .pem that we provided to you. For instance if it was x45x.pem, save it as x45x.ppk

Great, now your key file is ready and we can start Putty. In Putty, enter your user name and the IP address in the format <user_name>@. Make sure that 22 is given as Port and that SSH is selected.

Next, on the menu on the left, expand the “SSH” panel by clicking on the + and select “Auth”. Then, select your new putty format key file (.ppk) with Browse. Do NOT click on Open yet.

To make sure that you will not get logged out if you are inactive for a while, you can configure PuTTY to automatically send ‘keepalive’ signals at regular intervals to keep the session active. Click on Connection, then insert 180 to send keepalive signals every 3 minutes. Do NOT click on Open yet.

To avoid having to change the settings each time you log in, you can save the PuTTY session. Click on Session to get back to the basic options window. Once you are happy with all settings, give the session a name and click Save. Now you are ready to start the session with “Open”. The first time PuTTY will display a security alert dialog box that asks whether you trust the host you are connecting to. Click yes.

When you log in the next time, you can just click on the saved session and click load. If the IP address changed in the mean time (e.g. because we stopped the Amazon instance over night), you will need to replace the IP address by the new one. I would then recommend to Save the changed settings. Then simply click Open to start the session.

If the IP address did not change and you just want to login again, you can also right-click on the putty symbol in the taskbar (provided that you have pinned it to the taskbar) and select the session.

You can also tranfer files back and forth with Filezilla, a handy software to move files from a remote server such as the Amazon cloud or a cluster of your university.

Open Filezilla and choose Edit -> Settings. Next, choose SFTP and Add the .pem key file as indicated below and click OK. Finally, enter the IP address and the user name and when you hit enter, it should connect you. Next time, you can use the Quickconnect dropdown menu, as long as the IP address has not changed in the meantime. Now you will see the file directory system (folders) on your local computer on the left and your folders on the amazon cloud on the right. You can now just drag and drop files from one side to the other.
