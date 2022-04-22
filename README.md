# AM5_Phys
This is the repository that contains AM5 physics!
## Table of contents
* [Running AM5](README.md#running-am5)
* [Guidelines for Contributing](README.md#guidelines-for-contributing)
* [Using am5_phys in existing xmls](README.md#using-am5_phys-in-existing-xmls)
## Running AM5
To run this code please refer to the instuctions in the [AM5xml repo](https://gitlab.gfdl.noaa.gov/m5/am5xml#quickstart)
## Guidelines for Contributing
To contribute to this code. Please refer to the [contribution guide](https://gitlab.gfdl.noaa.gov/fms/am5_phys/-/blob/main/CONTRIBUTING.md)
## Using am5_phys in existing xmls
There are a few modifications to an XML in order switch to the am5_phys code. 
An XML will have a section in the compile experiment that looks like the following:
```xml
    <component name="atmos_phys" requires="fms" paths="atmos_phys">
      <description domainName="" communityName="" communityVersion="$(RELEASE)" communityGrid=""/>
      <source versionControl="git" root="http://gitlab.gfdl.noaa.gov/fms">
        <codeBase version="$(RELEASE)"> atmos_phys.git </codeBase>
          <csh><![CDATA[
            ( cd atmos_phys  && git checkout $(ATMOS_GIT_TAG) )
           ]]>
          </csh>
      </source>
      <compile>
        <cppDefs>$(F2003_FLAGS) -DCLUBB</cppDefs>
      </compile>
    </component>
```
1. The **paths** has been changed to "am5_phys".
2. The **codeBase** was changed from **atmos_phys.git** to **am5_phys.git**
3. The `csh` block was updated to only `cd` to am5_phys
