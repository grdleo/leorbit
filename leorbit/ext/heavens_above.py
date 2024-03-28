from typing import Optional
from bs4 import BeautifulSoup
from requests import get as get_url

def get_standard_magnitude(catnr: int) -> Optional[float]:
    url = f'https://www.heavens-above.com/SatInfo.aspx?satid={catnr}'

    try:
        res = get_url(url)
    except:
        raise Exception(f"Failed to connect to Heavens Above to retrieve standard magnitude for satellite")
    
    soup = BeautifulSoup(res.text, "html.parser")
    if "error" in soup.title.text.lower():
        raise ValueError(f"No satellite with CATNR {catnr}")
    
    bright = soup.find(id="ctl00_cph1_lblBrightness")
    
    if bright is None:
        return None # no brightness information
    brightness_txt: str
    try:
        brightness_txt = bright.contents[0].contents[1].text
    except:
        raise Exception("Heavens Above HTML might have changed")
    
    for s in brightness_txt.split(" "):
        try:
            return float(s)
        except:
            pass
    
    raise Exception("Heavens Above HTML might have changed")