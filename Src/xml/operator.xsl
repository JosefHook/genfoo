<?xml version="1.0" encoding="UTF-8"?>
<?modxslt-stylesheet type="text/xsl" media="fuffa, screen and $GET[stylesheet]" href="./%24GET%5Bstylesheet%5D" alternate="no" title="Translation using provided stylesheet" charset="ISO-8859-1" ?>
<?modxslt-stylesheet type="text/xsl" media="screen" alternate="no" title="Show raw source of the XML file" charset="ISO-8859-1" ?>

<!--
Stylesheet for generating GenFoo operator.hpp based on an xml description.

/author Josef Hook and Thomas Johnson
-->
<xsl:stylesheet 
   xmlns:yaslt="http://www.mod-xslt2.com/ns/1.0"
   xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"
   xmlns:xs="http://www.w3.org/2001/XMLSchema" 
   xmlns:fn="http://www.w3.org/2005/02/xpath-functions"
   xmlns:exsl="http://exslt.org/common"
   xmlns:str="http://exslt.org/strings"
   xmlns:func="http://exslt.org/functions"
   xmlns:my="http://localhost.localdomain/localns"
   exclude-result-prefixes="my"
   extension-element-prefixes="yaslt exsl func str">

<xsl:output method="text" version="1.0" encoding="UTF-8" indent="yes"/>

<!--
====================================================
Main template
====================================================
-->
    <xsl:template match="/Operator">

	    <xsl:variable name="className" select="'xmlDefinedOperator'"/>
	    <xsl:variable name="Dimensions" select="normalize-space(Dimensions)"/>

	    <xsl:text>
//
// Autogenerated from xml file operator.xml using the schema operator.xsd and the stylesheet operator.xsl.
//

#include &lt;cmath&gt;
#include &lt;libxml/tree.h&gt;
#include &lt;Operator.hpp&gt;
#include &lt;iostream&gt;
#include &lt;vector&gt;

</xsl:text>
<xsl:for-each select="Include">
  <xsl:text>
#include "</xsl:text>
  <xsl:value-of select="name"/>
  <xsl:text>"</xsl:text>
</xsl:for-each>
<xsl:text>

#define two_d_indexing(i,j,N) i*N+j

using namespace std;

// ==============================================================================
/**
</xsl:text><xsl:value-of select="Documentation"/><xsl:text>
*/
// ==============================================================================

class </xsl:text><xsl:value-of select="$className"/><xsl:text> : public Operator 
{

  public:
     </xsl:text><xsl:value-of select="$className"/><xsl:text>();
     void evalDrift(double* val, const double* x) const;
     void evalDiffusion(double* val, const double* x) const;
     void evalIC(double* val, const double* x) const;
     void evalSource(double* val, const double* x) const;

     void evalDriftEnumTerms(double* val, const double* x, int termIndex) const;
     void evalDiffusionEnumTerms(double* val, const double* x, int termIndex) const;
     void evalSourceEnumTerms(double* val, const double* x, int termIndex) const;

     int numberDriftTerms() const;
     int numberDiffusionTerms() const;
     int numberSourceTerms() const;

  private:
</xsl:text>
  <xsl:call-template name="declareVariables">
	  <xsl:with-param name="typeStr" select="type"/>
  </xsl:call-template>
  <xsl:call-template name="declareFunctions"/>
<xsl:text>
};
// End of Class

// ==============================================================================
// IMPLEMENTATIONS
// ==============================================================================
//
// The constructor
</xsl:text><xsl:value-of select="$className"/><xsl:text>::</xsl:text><xsl:value-of select="$className"/><xsl:text>() {
   _dim = </xsl:text><xsl:value-of select="$Dimensions"/><xsl:text>;
   _info = "</xsl:text><xsl:value-of select="normalize-space(Documentation)"/><xsl:text>";


</xsl:text>
<xsl:call-template name="implementVectorVariables">
  <xsl:with-param name="className" select="$className"/>
</xsl:call-template>
<xsl:text>
}

</xsl:text>
<!-- =================================================== -->
<xsl:for-each select="Drift">
<xsl:text>
/** </xsl:text><xsl:value-of select="Documentation"/>
<xsl:text>
*/
void </xsl:text><xsl:value-of select="$className"/><xsl:text>::evalDrift(double* val , const double* x) const {
</xsl:text>
    <xsl:value-of select="code"/>
    <xsl:text>
	    for (int j=0; j&lt; _dim   ; j++ ) {
               val[j] = 0.0;
	    }</xsl:text>
	    <xsl:for-each select="component"><xsl:text>
	    val[</xsl:text><xsl:value-of select="normalize-space(@index)"/><xsl:text>] = </xsl:text>
	    <xsl:value-of select="normalize-space(value)"/><xsl:text>;</xsl:text>
    </xsl:for-each>
    <xsl:text>
     }
</xsl:text>
</xsl:for-each>
<xsl:text>

</xsl:text>
<!-- =================================================== -->
<xsl:for-each select="Diffusion">
<xsl:text>
/** </xsl:text><xsl:value-of select="Documentation"/>
<xsl:text>
*/
void </xsl:text><xsl:value-of select="$className"/><xsl:text>::evalDiffusion(double* val , const double* x) const {
</xsl:text>
    <xsl:value-of select="code"/>
    <xsl:text>
	    for (int j=0; j&lt; _dim*_dim   ; j++ ) {
               val[j] = 0.0;
	    }</xsl:text>

<!-- std::cout &lt;&lt;"New:"&lt;&lt; two_d_indexing(0,0,2)&lt;&lt;  two_d_indexing(0,1,2)&lt;&lt; two_d_indexing(1,0,2)&lt;&lt; two_d_indexing(1,1,2)&lt;&lt;"; ";  -->

    <xsl:for-each select="component"><xsl:text>
	    val[two_d_indexing(</xsl:text><xsl:value-of select="normalize-space(@indexColumn)"/><xsl:text>,</xsl:text>
	    <xsl:value-of select="normalize-space(@indexRow)"/><xsl:text>,</xsl:text><xsl:value-of select="$Dimensions"/><xsl:text>)] = </xsl:text>
	    <xsl:value-of select="normalize-space(value)"/><xsl:text>;</xsl:text>
    </xsl:for-each>
<xsl:text>
}
</xsl:text>
</xsl:for-each>
<xsl:text>

</xsl:text>
<!-- =================================================== -->
<xsl:for-each select="Source">
<xsl:text>
/** </xsl:text><xsl:value-of select="Documentation"/>
<xsl:text>
*/
void </xsl:text><xsl:value-of select="$className"/><xsl:text>::evalSource(double* val , const double* x) const {
            </xsl:text><xsl:value-of select="code"/><xsl:text>
            val[0] = </xsl:text><xsl:value-of select="normalize-space(value)"/><xsl:text>;</xsl:text>
<xsl:text>
}
</xsl:text>
</xsl:for-each>
<xsl:text>

</xsl:text>
<!-- =================================================== -->
<xsl:for-each select="InitialCondition">
<xsl:text>
/** </xsl:text><xsl:value-of select="Documentation"/>
<xsl:text>
*/
void </xsl:text><xsl:value-of select="$className"/><xsl:text>::evalIC(double* val , const double* x) const {
            </xsl:text><xsl:value-of select="code"/><xsl:text>
            val[0] = </xsl:text><xsl:value-of select="normalize-space(value)"/><xsl:text>;</xsl:text>
<xsl:text>
}
</xsl:text>
</xsl:for-each>
<xsl:text>

// ==============================================================================
// Private variables
// ==============================================================================
</xsl:text>
<xsl:call-template name="implementScalarVariables">
  <xsl:with-param name="className" select="$className"/>
</xsl:call-template>
<xsl:text>
// ==============================================================================
// Private methods
// ==============================================================================
</xsl:text>
<xsl:call-template name="implementFunctions">
  <xsl:with-param name="className" select="$className"/>
</xsl:call-template>
<xsl:text>

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
Operator *load()  { return new </xsl:text><xsl:value-of select="$className"/><xsl:text>;  }
}

</xsl:text>
</xsl:template>
<!--
====================================================
END: Main template
====================================================
-->

<!--
====================================================
VARIABLES: Declare local/private variables
====================================================
-->
<xsl:template name="declareVariables">
  <xsl:for-each select="Private/variable">
    <xsl:choose>
      <xsl:when test="@isVector">
	<!-- VECTORS -->
	<xsl:text>     vector&lt;</xsl:text>
	<xsl:call-template name="cppDeclareType">
	  <xsl:with-param name="type" select="@type"/>
	  <xsl:with-param name="noSpaceAfterTypeDecleration" select="1"/>
	</xsl:call-template>
	<xsl:text>&gt; </xsl:text>
	<xsl:value-of select="normalize-space(name)"/><xsl:text>;
</xsl:text>
      </xsl:when>
      <xsl:otherwise>
	<!-- SCALARS -->
	<xsl:text>     </xsl:text>
	<xsl:if test="@static">
	  <xsl:text>static </xsl:text>
	</xsl:if>
	<xsl:if test="@const">
	  <xsl:text>const </xsl:text>
	</xsl:if>
	<xsl:call-template name="cppDeclareType">
	  <xsl:with-param name="type" select="@type"/>
	</xsl:call-template>
	<xsl:value-of select="normalize-space(name)"/><xsl:text>;
</xsl:text>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:for-each>
</xsl:template>


<!--
====================================================
VARIABLES: Implement scalar private variables
====================================================
-->
<xsl:template name="implementScalarVariables">
  <xsl:param name="className"/>
  <xsl:for-each select="Private/variable">
    <xsl:choose>
      <xsl:when test="not(@isVector)">
	<xsl:text>
/** </xsl:text><xsl:value-of select="documentation"/>
        <xsl:if test="unit">
	  <xsl:text>
Unit: </xsl:text><xsl:value-of select="normalize-space(unit)"/>
        </xsl:if>
	<xsl:text>
*/
</xsl:text>
	<xsl:if test="@const and value">
	  <xsl:text>const </xsl:text>
          <xsl:call-template name="cppDeclareType">
	    <xsl:with-param name="type" select="@type"/>
	  </xsl:call-template>
	  <xsl:value-of select="$className"/><xsl:text>::</xsl:text>
	  <xsl:value-of select="normalize-space(name)"/><xsl:text>=</xsl:text>
	  <xsl:value-of select="normalize-space(value)"/><xsl:text>;
</xsl:text>
	</xsl:if>
      </xsl:when>
    </xsl:choose>
  </xsl:for-each>
</xsl:template>


<!--
====================================================
VARIABLES: Implement vector private variables
====================================================
-->
<xsl:template name="implementVectorVariables">
  <xsl:param name="className"/>
  <xsl:for-each select="Private/variable">
    <xsl:if test="@isVector">
      <xsl:text>
  /** </xsl:text><xsl:value-of select="documentation"/>
        <xsl:if test="unit">
	  <xsl:text>
  Unit: </xsl:text><xsl:value-of select="normalize-space(unit)"/>
        </xsl:if>
        <xsl:text>
  */
  </xsl:text>
	<xsl:value-of select="normalize-space(name)"/><xsl:text>.clear();
</xsl:text>
        <xsl:call-template name="pushNextCommaSeparatedValue">
	  <xsl:with-param name="valueList" select="value"/>
	</xsl:call-template>
      </xsl:if>
  </xsl:for-each>
</xsl:template>

<!--
====================================================
GENERIC: Push next value from comma-separated list into vector
====================================================
-->
<xsl:template name="pushNextCommaSeparatedValue">
  <xsl:param name="valueList"/>
  <xsl:text>  </xsl:text>
  <xsl:value-of select="normalize-space(name)"/>
  <xsl:text>.push_back(</xsl:text>
  <xsl:choose>
    <xsl:when test="contains($valueList,',')">
      <xsl:value-of select="substring-before($valueList,',')"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="$valueList"/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>);
</xsl:text>
  <xsl:if test="string-length(substring-before($valueList,','))>0">
  <xsl:call-template name="pushNextCommaSeparatedValue">
    <xsl:with-param name="valueList" select="substring-after($valueList,',')"/>
  </xsl:call-template>
  </xsl:if>
</xsl:template>


<!--
====================================================
FUNCTIONS: Declare functions
====================================================
-->
<xsl:template name="declareFunctions">
  <xsl:for-each select="Private/function">
    <xsl:text>     </xsl:text>
    <xsl:call-template name="cppDeclareType">
      <xsl:with-param name="type" select="@type"/>
    </xsl:call-template>
    <xsl:value-of select="normalize-space(name)"/>
    <xsl:text>(</xsl:text>
    <xsl:value-of select="input"/>
    <xsl:text>) const;
</xsl:text>
  </xsl:for-each>
</xsl:template>



<!--
====================================================
FUNCTIONS: Implement functions
====================================================
-->
<xsl:template name="implementFunctions">
  <xsl:param name="className"/>

  <xsl:for-each select="Private/function">
  <xsl:call-template name="implementOneFunction">
    <xsl:with-param name="className" select="$className"/>
    <xsl:with-param name="type" select="@type"/>
    <xsl:with-param name="name" select="name"/>
    <xsl:with-param name="input" select="input"/>
    <xsl:with-param name="code" select="code"/>
    <xsl:with-param name="documentation" select="documentation"/>
    <xsl:with-param name="unit" select="unit"/>
  </xsl:call-template>
  </xsl:for-each>
</xsl:template>

<!--
====================================================
FUNCTIONS: Implement functions
====================================================
-->
<xsl:template name="implementOneFunction">
  <xsl:param name="className"/>
  <xsl:param name="type"/>
  <xsl:param name="name"/>
  <xsl:param name="input"/>
  <xsl:param name="code"/>
  <xsl:param name="code1"/>
  <xsl:param name="documentation"/>
  <xsl:param name="unit"/>

    <xsl:text>
/** </xsl:text><xsl:value-of select="$documentation"/>
    <xsl:if test="$unit">
      <xsl:text>
Unit: </xsl:text><xsl:value-of select="normalize-space($unit)"/>
    </xsl:if>
    <xsl:text>
*/
</xsl:text>
    <xsl:call-template name="cppDeclareType">
      <xsl:with-param name="type" select="$type"/>
    </xsl:call-template>
    <xsl:value-of select="$className"/>
    <xsl:text>::</xsl:text>
    <xsl:value-of select="normalize-space($name)"/>
    <xsl:text>(</xsl:text>
    <xsl:value-of select="$input"/>
    <xsl:text>) const {
</xsl:text>
    <xsl:value-of select="$code"/>
    <xsl:text>
</xsl:text>
    <xsl:value-of select="$code1"/>
    <xsl:text>
}
</xsl:text>
</xsl:template>



<!--
====================================================
GENERIC: Declare a variable type.
Translate input $typeStr into an intrinsic type.
====================================================
-->
<xsl:template name="cppDeclareType">
  <xsl:param name="type"/>
  <xsl:param name="noSpaceAfterTypeDecleration"/>

  <xsl:choose>
    <xsl:when test="$type='float'">
      <xsl:text>double</xsl:text>
    </xsl:when>
    <xsl:when test="$type='double'">
      <xsl:text>double</xsl:text>
    </xsl:when>
    <xsl:when test="$type='decimal'">
      <xsl:text>double</xsl:text>
    </xsl:when>
    <xsl:when test="$type='real'">
      <xsl:text>double</xsl:text>
    </xsl:when>
    <xsl:when test="$type='int'">
      <xsl:text>int</xsl:text>
    </xsl:when>
    <xsl:when test="$type='integer'">
      <xsl:text>int</xsl:text>
    </xsl:when>
    <xsl:when test="$type='boolean'">
      <xsl:text>boolean</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="normalize-space($type)"/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:if test="not($noSpaceAfterTypeDecleration)">
    <xsl:text> </xsl:text>
  </xsl:if>
</xsl:template>

</xsl:stylesheet>
